#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <fstream>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/search/search.hpp>

#include <cereal/archives/binary.hpp>

// -------- unified "dna5_vector" container type (portable) --------
using dna5_vec = std::vector<seqan3::dna5>;
using dna5_collection = std::vector<dna5_vec>;
// ----------------------------------------------------------------

// Convert "whatever SeqAn3 returns" into two numbers (a,b).
template <typename T>
static inline std::pair<size_t, size_t> to_two_numbers(T const &v)
{
    if constexpr (requires { v.first; v.second; }) // pair-like
    {
        return {static_cast<size_t>(v.first), static_cast<size_t>(v.second)};
    }
    else if constexpr (requires { std::get<0>(v); std::get<1>(v); }) // tuple-like
    {
        return {static_cast<size_t>(std::get<0>(v)), static_cast<size_t>(std::get<1>(v))};
    }
    else if constexpr (std::is_integral_v<std::remove_reference_t<T>>) // single number
    {
        return {static_cast<size_t>(v), 0u}; // note: interpret as "pos only" later
    }
    else
    {
        static_assert(!sizeof(T), "Unsupported reference position type returned by SeqAn3 search result.");
    }
}

template <typename SearchResult>
static inline auto get_ref_begin_any(SearchResult const &res)
{
    if constexpr (requires { res.reference_begin_position(); })
        return res.reference_begin_position();
    else if constexpr (requires { res.reference_begin(); })
        return res.reference_begin();
    else
        static_assert(!sizeof(SearchResult), "No reference_begin_position()/reference_begin() in result.");
}

// Robustly decode (text_id, pos) for collection + single-text across SeqAn3 versions.
template <typename SearchResult>
static inline bool decode_text_id_and_pos(SearchResult const &res,
                                         dna5_collection const &reference_texts,
                                         size_t &text_id,
                                         size_t &pos)
{
    size_t ntexts = reference_texts.size();
    if (ntexts == 0) return false;

    // 1) If SeqAn3 provides reference_id(), use it. This is the most reliable for collections.
    if constexpr (requires { res.reference_id(); })
    {
        text_id = static_cast<size_t>(res.reference_id());
        if (text_id >= ntexts) return false;

        auto const &v = get_ref_begin_any(res);

        // v might be an integer "pos", or a pair/tuple containing pos.
        if constexpr (std::is_integral_v<std::remove_reference_t<decltype(v)>>)
        {
            pos = static_cast<size_t>(v);
            return pos < reference_texts[text_id].size();
        }
        else
        {
            auto [a, b] = to_two_numbers(v);
            // choose which one is the pos by checking range in this specific text
            bool a_ok = (a < reference_texts[text_id].size());
            bool b_ok = (b < reference_texts[text_id].size());

            if (a_ok && !b_ok) { pos = a; return true; }
            if (b_ok && !a_ok) { pos = b; return true; }
            if (a_ok && b_ok)  { pos = std::min(a, b); return true; } // both small -> take smaller
            return false;
        }
    }
    else
    {
        // 2) Fallback: no reference_id(). Then reference_begin must encode (text_id,pos) somehow.
        auto const &v = get_ref_begin_any(res);

        if constexpr (std::is_integral_v<std::remove_reference_t<decltype(v)>>)
        {
            // If we only get a pos and no reference_id, we can only assume single-text.
            text_id = 0;
            pos = static_cast<size_t>(v);
            return (pos < reference_texts[0].size());
        }
        else
        {
            auto [a, b] = to_two_numbers(v);
            bool ab_ok = (a < ntexts) && (b < reference_texts[a].size());
            bool ba_ok = (b < ntexts) && (a < reference_texts[b].size());

            if (ab_ok && !ba_ok) { text_id = a; pos = b; return true; }
            if (ba_ok && !ab_ok) { text_id = b; pos = a; return true; }
            if (ab_ok && ba_ok)
            {
                if (a <= b) { text_id = a; pos = b; }
                else        { text_id = b; pos = a; }
                return true;
            }
            return false;
        }
    }
}

// Hamming distance <= k between query q and reference ref at offset start (no copies).
static inline bool hamming_ref_leq_k(dna5_vec const &q, dna5_vec const &ref, int start, int k)
{
    int mism = 0;
    int m = static_cast<int>(q.size());
    for (int i = 0; i < m; ++i)
    {
        mism += (seqan3::to_char(q[static_cast<size_t>(i)]) !=
                 seqan3::to_char(ref[static_cast<size_t>(start + i)]));
        if (mism > k)
            return false;
    }
    return true;
}

int main(int argc, char const *const *argv)
{
    seqan3::argument_parser parser{"fmindex_pigeon_search", argc, argv, seqan3::update_notifications::off};
    parser.info.author = "SeqAn-Team";
    parser.info.version = "1.0.0";

    std::filesystem::path index_path{};
    parser.add_option(index_path, '\0', "index", "Path to the serialized FM-index (.bin).");

    std::filesystem::path reference_file{};
    parser.add_option(reference_file, '\0', "reference", "Path to the reference FASTA/FASTQ file (for verification).");

    std::filesystem::path query_file{};
    parser.add_option(query_file, '\0', "query", "Path to the query FASTA/FASTQ file.");

    size_t number_of_queries{100};
    parser.add_option(number_of_queries, '\0', "query_ct",
                      "Number of queries; if not enough queries, they will be duplicated.");

    uint8_t number_of_errors{0};
    parser.add_option(number_of_errors, '\0', "errors",
                      "Maximum allowed Hamming errors (mismatches / substitutions only).");

    try { parser.parse(); }
    catch (seqan3::argument_parser_error const &ext)
    {
        seqan3::debug_stream << "Parsing error: " << ext.what() << "\n";
        return EXIT_FAILURE;
    }

    // Load reference texts
    seqan3::sequence_file_input reference_stream{reference_file};
    dna5_collection reference_texts{};
    for (auto &rec : reference_stream)
    {
        auto const &seq = rec.sequence();
        reference_texts.emplace_back(seq.begin(), seq.end());
    }
    if (reference_texts.empty())
    {
        seqan3::debug_stream << "Error: reference file contains no sequences.\n";
        return EXIT_FAILURE;
    }

    // Load queries
    seqan3::sequence_file_input query_stream{query_file};
    dna5_collection queries{};
    for (auto &rec : query_stream)
    {
        auto const &seq = rec.sequence();
        queries.emplace_back(seq.begin(), seq.end());
    }
    if (queries.empty() && number_of_queries > 0)
    {
        seqan3::debug_stream << "Error: query file contains no sequences.\n";
        return EXIT_FAILURE;
    }

    if (number_of_queries == 0) queries.clear();
    else
    {
        while (queries.size() < number_of_queries)
        {
            size_t old_count = queries.size();
            queries.resize(2 * old_count);
            std::copy_n(queries.begin(), old_count, queries.begin() + old_count);
        }
        queries.resize(number_of_queries);
    }

    // Load index
    using Index = decltype(seqan3::fm_index{dna5_collection{}});
    Index index{};
    try
    {
        std::ifstream is{index_path, std::ios::binary};
        if (!is)
        {
            seqan3::debug_stream << "Error: cannot open index file: " << index_path << "\n";
            return EXIT_FAILURE;
        }
        cereal::BinaryInputArchive iarchive{is};
        iarchive(index);
    }
    catch (std::exception const &e)
    {
        seqan3::debug_stream << "Error while loading index: " << e.what() << "\n";
        return EXIT_FAILURE;
    }

    int k = static_cast<int>(number_of_errors);

    // Exact search for seeds (0 error)
    seqan3::configuration const cfg_exact =
        seqan3::search_cfg::max_error_substitution{seqan3::search_cfg::error_count{0}} |
        seqan3::search_cfg::max_error_insertion{seqan3::search_cfg::error_count{0}} |
        seqan3::search_cfg::max_error_deletion{seqan3::search_cfg::error_count{0}};

    size_t hits = 0;
    auto t0 = std::chrono::high_resolution_clock::now();

    for (auto const &q : queries)
    {
        int m = static_cast<int>(q.size());
        if (m == 0) continue;

        int parts = k + 1;
        if (parts < 1) parts = 1;
        if (parts > m) parts = m;

        std::vector<int> cut(parts + 1, 0);
        for (int i = 0; i <= parts; ++i)
            cut[i] = static_cast<int>((static_cast<long long>(i) * m) / parts);

        std::vector<std::pair<size_t, int>> cand;
        cand.reserve(256);

        for (int p = 0; p < parts; ++p)
        {
            int qs = cut[p];
            int qe = cut[p + 1];
            if (qe <= qs) continue;

            dna5_vec piece(q.begin() + qs, q.begin() + qe);

            for (auto &&res : seqan3::search(piece, index, cfg_exact))
            {
                size_t text_id = 0, pos = 0;
                if (!decode_text_id_and_pos(res, reference_texts, text_id, pos))
                    continue;

                int start = static_cast<int>(pos) - qs;
                cand.push_back({text_id, start});
            }
        }

        if (cand.empty()) continue;

        std::sort(cand.begin(), cand.end());
        cand.erase(std::unique(cand.begin(), cand.end()), cand.end());

        for (auto const &[text_id, start] : cand)
        {
            auto const &ref = reference_texts[text_id];
            if (start < 0) continue;
            if (start + m > static_cast<int>(ref.size())) continue;

            if (hamming_ref_leq_k(q, ref, start, k))
                ++hits;
        }
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    double search_s = std::chrono::duration<double>(t1 - t0).count();

    seqan3::debug_stream << "Search time: " << search_s << " seconds.\n";
    seqan3::debug_stream << "queries=" << queries.size()
                         << " errors=" << static_cast<int>(number_of_errors)
                         << " hits=" << hits << "\n";
    return EXIT_SUCCESS;
}

