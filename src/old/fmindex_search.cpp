#include <filesystem>
#include <vector>
#include <algorithm>
#include <cstdint>
#include <chrono>
#include <thread>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/concept.hpp> // assign_char_to
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/search/search.hpp>

using dna5_vec = std::vector<seqan3::dna5>;

int main(int argc, char const *const *argv)
{
    seqan3::argument_parser parser{"fmindex_search", argc, argv, seqan3::update_notifications::off};
    parser.info.author = "SeqAn-Team";
    parser.info.version = "1.0.0";

    std::filesystem::path reference_file{};
    parser.add_option(reference_file, '\0', "reference", "Path to the reference FASTA/FASTQ file.");

    std::filesystem::path query_file{};
    parser.add_option(query_file, '\0', "query", "Path to the query FASTA/FASTQ file.");

    size_t number_of_queries{100};
    parser.add_option(number_of_queries, '\0', "query_ct",
                      "Number of queries; if not enough, queries will be duplicated logically (no extra memory).");

    uint8_t number_of_errors{0};
    parser.add_option(number_of_errors, '\0', "errors", "Allowed substitutions (Hamming distance).");

    size_t threads{0};
    parser.add_option(threads, '\0', "threads", "Number of worker threads (0 = use hardware_concurrency).");

    size_t guard{50};
    parser.add_option(guard, '\0', "guard", "Number of 'N' separators between reference records (default 50).");

    try
    {
        parser.parse();
    }
    catch (seqan3::argument_parser_error const &e)
    {
        seqan3::debug_stream << "Parsing error: " << e.what() << "\n";
        return EXIT_FAILURE;
    }

    // -----------------------------
    // Load reference into ONE dna5_vec
    // -----------------------------
    seqan3::sequence_file_input reference_stream{reference_file};

    dna5_vec reference_concat{};
    reference_concat.reserve(1'000'000);

    seqan3::dna5 sep{};
    seqan3::assign_char_to('N', sep);

    size_t ref_seq_count = 0;
    size_t ref_total_len = 0;

    for (auto &&rec : reference_stream)
    {
        auto const &seq = rec.sequence();
        if (!seq.empty())
        {
            reference_concat.insert(reference_concat.end(), seq.begin(), seq.end());
            for (size_t i = 0; i < guard; ++i)
                reference_concat.push_back(sep);
            ref_total_len += seq.size();
        }
        ++ref_seq_count;
    }

    if (ref_seq_count == 0 || reference_concat.empty())
    {
        seqan3::debug_stream << "Error: reference file contains no sequences.\n";
        return EXIT_FAILURE;
    }

    // remove trailing guard
    for (size_t i = 0; i < guard && !reference_concat.empty(); ++i)
        reference_concat.pop_back();

    // -----------------------------
    // Load base queries (NO DUP COPIES)
    // -----------------------------
    seqan3::sequence_file_input query_stream{query_file};
    std::vector<dna5_vec> base_queries{};
    base_queries.reserve(256);

    for (auto &&rec : query_stream)
        base_queries.push_back(rec.sequence());

    if (base_queries.empty())
    {
        seqan3::debug_stream << "Error: query file contains no sequences.\n";
        return EXIT_FAILURE;
    }

    if (number_of_queries == 0)
    {
        seqan3::debug_stream << "queries=0 errors=" << static_cast<uint32_t>(number_of_errors)
                             << " threads=0 hits=0\n";
        return 0;
    }

    // threads
    if (threads == 0)
    {
        threads = std::thread::hardware_concurrency();
        if (threads == 0)
            threads = 1;
    }
    threads = std::min(threads, number_of_queries);
    if (threads == 0)
        threads = 1;

    // -----------------------------
    // Build index (MOVE in)
    // -----------------------------
    auto t_index0 = std::chrono::steady_clock::now();
    seqan3::fm_index index{std::move(reference_concat)};
    auto t_index1 = std::chrono::steady_clock::now();

    seqan3::debug_stream << "Index Construction time: "
                         << std::chrono::duration<double>(t_index1 - t_index0).count()
                         << " seconds.\n";

    // force release moved-from vector capacity
    dna5_vec().swap(reference_concat);

    // -----------------------------
    // COUNT-ONLY search cfg (NO LOCATE OUTPUT)
    // -----------------------------
    seqan3::configuration const cfg =
        seqan3::search_cfg::hit_all{} |
        seqan3::search_cfg::max_error_substitution{seqan3::search_cfg::error_count{number_of_errors}} |
        seqan3::search_cfg::max_error_insertion{seqan3::search_cfg::error_count{0}} |
        seqan3::search_cfg::max_error_deletion{seqan3::search_cfg::error_count{0}};

    // -----------------------------
    // Parallel search
    // -----------------------------
    auto t0 = std::chrono::steady_clock::now();

    std::vector<std::thread> pool;
    pool.reserve(threads);
    std::vector<size_t> local_hits(threads, 0);

    auto worker = [&](size_t tid, size_t begin, size_t end)
    {
        size_t occ = 0;
        size_t const m = base_queries.size();

        for (size_t i = begin; i < end; ++i)
        {
            auto const &q = base_queries[i % m];
            for (auto &&res : seqan3::search(q, index, cfg))
            {
                (void)res;
                ++occ;
            }
        }
        local_hits[tid] = occ;
    };

    size_t n = number_of_queries;
    size_t block = (n + threads - 1) / threads;

    for (size_t t = 0; t < threads; ++t)
    {
        size_t begin = t * block;
        size_t end = std::min(n, begin + block);
        if (begin < end)
            pool.emplace_back(worker, t, begin, end);
    }

    for (auto &th : pool)
        th.join();

    size_t total_hits = 0;
    for (size_t h : local_hits)
        total_hits += h;

    auto t1 = std::chrono::steady_clock::now();
    seqan3::debug_stream << "Search time: "
                         << std::chrono::duration<double>(t1 - t0).count()
                         << " seconds.\n";

    seqan3::debug_stream << "queries=" << number_of_queries
                         << " base_queries=" << base_queries.size()
                         << " errors=" << static_cast<uint32_t>(number_of_errors)
                         << " threads=" << threads
                         << " hits=" << total_hits << "\n";

    return 0;
}
