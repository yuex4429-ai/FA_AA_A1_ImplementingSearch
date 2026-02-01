#include <sstream>
#include <filesystem>
#include <vector>
#include <algorithm>
#include <cstdint>
#include <unordered_set>
#include <chrono>
#include <thread>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/search/search.hpp>

// ---------- helpers ----------
struct candidate_key
{
    uint32_t ref_id{};
    uint64_t start{};
    bool operator==(candidate_key const &o) const noexcept { return ref_id == o.ref_id && start == o.start; }
};

struct candidate_key_hash
{
    size_t operator()(candidate_key const &k) const noexcept
    {
        size_t h1 = std::hash<uint32_t>{}(k.ref_id);
        size_t h2 = std::hash<uint64_t>{}(k.start);
        return h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1 << 6) + (h1 >> 2));
    }
};

static inline bool verify_hamming(std::vector<seqan3::dna5> const &text,
                                  std::vector<seqan3::dna5> const &query,
                                  int64_t start_pos,
                                  uint8_t max_errors) noexcept
{
    if (start_pos < 0)
        return false;
    uint64_t start = static_cast<uint64_t>(start_pos);
    if (start + query.size() > text.size())
        return false;

    uint8_t errors = 0;
    for (size_t i = 0; i < query.size(); ++i)
    {
        if (text[start + i] != query[i])
            if (++errors > max_errors)
                return false;
    }
    return true;
}

int main(int argc, char const *const *argv)
{
    seqan3::argument_parser parser{"fmindex_pigeon_search", argc, argv, seqan3::update_notifications::off};

    parser.info.author = "SeqAn-Team";
    parser.info.version = "1.0.0";

    std::filesystem::path reference_file{};
    parser.add_option(reference_file, '\0', "reference", "Path to the reference FASTA/FASTQ file.");

    std::filesystem::path query_file{};
    parser.add_option(query_file, '\0', "query", "Path to the query FASTA/FASTQ file.");

    size_t number_of_queries{100};
    parser.add_option(number_of_queries, '\0', "query_ct", "Number of queries; if not enough, queries will be duplicated.");

    uint8_t number_of_errors{0};
    parser.add_option(number_of_errors, '\0', "errors", "Allowed Hamming distance (substitutions only).");

    size_t threads{0};
    parser.add_option(threads, '\0', "threads", "Number of worker threads (0 = use hardware_concurrency).");

    try
    {
        parser.parse();
    }
    catch (seqan3::argument_parser_error const &e)
    {
        seqan3::debug_stream << "Parsing error: " << e.what() << "\n";
        return EXIT_FAILURE;
    }

    // load reference + queries
    seqan3::sequence_file_input reference_stream{reference_file};
    seqan3::sequence_file_input query_stream{query_file};

    std::vector<std::vector<seqan3::dna5>> reference{};
    for (auto &rec : reference_stream)
        reference.push_back(rec.sequence());

    std::vector<std::vector<seqan3::dna5>> queries{};
    for (auto &rec : query_stream)
        queries.push_back(rec.sequence());

    if (reference.empty() || queries.empty())
    {
        seqan3::debug_stream << "Error: empty reference or query file.\n";
        return EXIT_FAILURE;
    }

    while (queries.size() < number_of_queries)
    {
        size_t old = queries.size();
        queries.resize(2 * old);
        std::copy_n(queries.begin(), old, queries.begin() + old);
    }
    queries.resize(number_of_queries);

    if (threads == 0)
    {
        threads = std::thread::hardware_concurrency();
        if (threads == 0)
            threads = 1;
    }
    if (threads > queries.size())
        threads = queries.size();
    if (threads == 0)
        threads = 1;

    // build index timing
    auto t_index0 = std::chrono::steady_clock::now();
    seqan3::fm_index index{reference};
    auto t_index1 = std::chrono::steady_clock::now();
    double index_seconds = std::chrono::duration<double>(t_index1 - t_index0).count();
    seqan3::debug_stream << "Index Construction time: " << index_seconds << " seconds.\n";

    seqan3::configuration const exact_cfg =
        seqan3::search_cfg::hit_all{} |
        seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{0}} |
        seqan3::search_cfg::output_query_id{} |
        seqan3::search_cfg::output_reference_id{} |
        seqan3::search_cfg::output_reference_begin_position{};

    auto t0 = std::chrono::steady_clock::now();

    std::vector<std::thread> pool;
    pool.reserve(threads);
    std::vector<size_t> local_hits(threads, 0);

    auto worker = [&](size_t tid, size_t begin, size_t end)
    {
        size_t hits = 0;

        for (size_t qi = begin; qi < end; ++qi)
        {
            auto const &q = queries[qi];
            if (q.empty())
                continue;

            size_t k = static_cast<size_t>(number_of_errors) + 1;
            if (k > q.size())
                k = q.size();

            std::vector<size_t> b(k + 1, 0);
            for (size_t i = 0; i <= k; ++i)
                b[i] = (i * q.size()) / k;

            std::unordered_set<candidate_key, candidate_key_hash> seen;
            seen.reserve(64);

            for (size_t part_id = 0; part_id < k; ++part_id)
            {
                size_t pb = b[part_id];
                size_t pe = b[part_id + 1];
                if (pe <= pb)
                    continue;

                std::vector<seqan3::dna5> part;
                part.reserve(pe - pb);
                part.insert(part.end(), q.begin() + pb, q.begin() + pe);

                for (auto &&res : seqan3::search(part, index, exact_cfg))
                {
                    uint32_t ref_id = static_cast<uint32_t>(res.reference_id());
                    int64_t hit_pos = static_cast<int64_t>(res.reference_begin_position());
                    int64_t start = hit_pos - static_cast<int64_t>(pb);
                    if (start < 0)
                        continue;

                    candidate_key key{ref_id, static_cast<uint64_t>(start)};
                    if (!seen.insert(key).second)
                        continue;

                    if (ref_id < reference.size() && verify_hamming(reference[ref_id], q, start, number_of_errors))
                        ++hits;
                }
            }
        }

        local_hits[tid] = hits;
    };

    size_t n = queries.size();
    size_t block = (n + threads - 1) / threads;

    for (size_t t = 0; t < threads; ++t)
    {
        size_t begin = t * block;
        size_t end = std::min(n, begin + block);
        if (begin >= end)
        {
            local_hits[t] = 0;
            continue;
        }
        pool.emplace_back(worker, t, begin, end);
    }
    for (auto &th : pool)
        th.join();

    size_t total_hits = 0;
    for (size_t h : local_hits)
        total_hits += h;

    auto t1 = std::chrono::steady_clock::now();
    double search_seconds = std::chrono::duration<double>(t1 - t0).count();
    seqan3::debug_stream << "Search time: " << search_seconds << " seconds.\n";

    seqan3::debug_stream << "queries=" << queries.size()
                         << " errors=" << static_cast<uint32_t>(number_of_errors)
                         << " threads=" << threads
                         << " verified_hits=" << total_hits << "\n";

    return 0;
}
