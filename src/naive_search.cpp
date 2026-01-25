#include <filesystem>
#include <vector>
#include <string>
#include <algorithm>
#include <cstdint>
#include <chrono>
#include <thread>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>

// dna5 -> std::string (ACGTN)
static inline std::string to_string_dna5(std::vector<seqan3::dna5> const &seq)
{
    std::string s;
    s.resize(seq.size());
    for (size_t i = 0; i < seq.size(); ++i)
        s[i] = seqan3::to_char(seq[i]);
    return s;
}

// Duplicate queries to N (doubling like your Python/C++)
static inline void duplicate_to_n_cppstyle(std::vector<std::string> &queries, size_t n)
{
    if (n == 0 || queries.empty())
    {
        queries.clear();
        return;
    }
    while (queries.size() < n)
    {
        size_t old = queries.size();
        queries.resize(2 * old);
        std::copy_n(queries.begin(), old, queries.begin() + old);
    }
    queries.resize(n);
}

// Split [0,n) into blocks, at least min(threads,n) blocks, prefer >= min_block
static inline std::vector<std::pair<size_t, size_t>>
chunk_ranges(size_t n, size_t threads, size_t min_block)
{
    std::vector<std::pair<size_t, size_t>> rs;
    if (n == 0)
        return rs;

    threads = std::max<size_t>(1, threads);
    min_block = std::max<size_t>(1, min_block);

    size_t min_blocks = std::min(threads, n);
    size_t blocks_by_min_block = (n + min_block - 1) / min_block;
    size_t blocks = std::max(min_blocks, blocks_by_min_block);
    blocks = std::min(blocks, n);

    size_t block_size = (n + blocks - 1) / blocks;

    rs.reserve(blocks);
    for (size_t t = 0; t < blocks; ++t)
    {
        size_t b = t * block_size;
        size_t e = std::min(n, b + block_size);
        if (b < e)
            rs.emplace_back(b, e);
    }
    return rs;
}

// Count overlaps using std::string::find (usually highly optimized)
static inline size_t count_overlaps_find(std::string const &text, std::string const &pat)
{
    if (pat.empty() || pat.size() > text.size())
        return 0;

    size_t cnt = 0;
    size_t pos = 0;
    while (true)
    {
        pos = text.find(pat, pos);
        if (pos == std::string::npos)
            break;
        ++cnt;
        ++pos; // overlap allowed
    }
    return cnt;
}

int main(int argc, char const *const *argv)
{
    seqan3::argument_parser parser{"naive_search", argc, argv, seqan3::update_notifications::off};
    parser.info.author = "SeqAn-Team";
    parser.info.version = "1.0.0";

    std::filesystem::path reference_file{};
    parser.add_option(reference_file, '\0', "reference", "Path to the reference FASTA/FASTQ file.");

    std::filesystem::path query_file{};
    parser.add_option(query_file, '\0', "query", "Path to the query FASTA/FASTQ file.");

    size_t number_of_queries{100};
    parser.add_option(number_of_queries, '\0', "query_ct",
                      "Number of queries; if not enough, queries will be duplicated.");

    uint8_t number_of_errors{0};
    parser.add_option(number_of_errors, '\0', "errors",
                      "Allowed substitutions. NOTE: naive_search supports exact match only; errors forced to 0.");

    size_t threads{0};
    parser.add_option(threads, '\0', "threads",
                      "Number of worker threads (0 = use hardware_concurrency).");

    size_t min_block{256};
    parser.add_option(min_block, '\0', "min_block",
                      "Minimum number of queries per block (controls granularity).");

    try
    {
        parser.parse();
    }
    catch (seqan3::argument_parser_error const &e)
    {
        seqan3::debug_stream << "Parsing error: " << e.what() << "\n";
        return EXIT_FAILURE;
    }

    // exact only
    number_of_errors = 0;

    // Load reference (multi-chr supported)
    seqan3::sequence_file_input reference_stream{reference_file};
    std::vector<std::string> reference;
    reference.reserve(64);
    for (auto &rec : reference_stream)
        reference.push_back(to_string_dna5(rec.sequence()));

    if (reference.empty())
    {
        seqan3::debug_stream << "Error: reference file contains no sequences.\n";
        return EXIT_FAILURE;
    }

    // Load queries
    seqan3::sequence_file_input query_stream{query_file};
    std::vector<std::string> queries;
    for (auto &rec : query_stream)
        queries.push_back(to_string_dna5(rec.sequence()));

    if (queries.empty())
    {
        seqan3::debug_stream << "Error: query file contains no sequences.\n";
        return EXIT_FAILURE;
    }

    duplicate_to_n_cppstyle(queries, number_of_queries);
    if (queries.empty())
    {
        seqan3::debug_stream << "Error: No queries after duplication.\n";
        return EXIT_FAILURE;
    }

    // threads
    if (threads == 0)
    {
        threads = std::thread::hardware_concurrency();
        if (threads == 0)
            threads = 1;
    }
    threads = std::min(threads, queries.size());
    if (threads == 0)
        threads = 1;

    auto ranges = chunk_ranges(queries.size(), threads, min_block);
    size_t used_threads = std::min(threads, ranges.size());
    if (used_threads == 0)
        used_threads = 1;

    // Search timing ONLY (no extra preprocessing inside)
    auto t0 = std::chrono::steady_clock::now();

    std::vector<std::thread> pool;
    pool.reserve(used_threads);
    std::vector<size_t> local_hits(used_threads, 0);

    auto worker = [&](size_t tid, size_t b, size_t e)
    {
        size_t hits = 0;
        for (size_t qi = b; qi < e; ++qi)
        {
            auto const &q = queries[qi];
            if (q.empty())
                continue;

            for (auto const &chr : reference)
                hits += count_overlaps_find(chr, q);
        }
        local_hits[tid] = hits;
    };

    for (size_t t = 0; t < used_threads; ++t)
    {
        auto [b, e] = ranges[t];
        pool.emplace_back(worker, t, b, e);
    }

    size_t main_hits = 0;
    for (size_t t = used_threads; t < ranges.size(); ++t)
    {
        auto [b, e] = ranges[t];
        for (size_t qi = b; qi < e; ++qi)
        {
            auto const &q = queries[qi];
            if (q.empty())
                continue;
            for (auto const &chr : reference)
                main_hits += count_overlaps_find(chr, q);
        }
    }

    for (auto &th : pool)
        th.join();

    size_t total_hits = main_hits;
    for (size_t t = 0; t < used_threads; ++t)
        total_hits += local_hits[t];

    auto t1 = std::chrono::steady_clock::now();
    double seconds = std::chrono::duration<double>(t1 - t0).count();

    seqan3::debug_stream << "Search time: " << seconds << " seconds.\n";
    seqan3::debug_stream << "queries=" << queries.size()
                         << " errors=" << static_cast<uint32_t>(number_of_errors)
                         << " threads=" << used_threads
                         << " hits=" << total_hits << "\n";

    return 0;
}
