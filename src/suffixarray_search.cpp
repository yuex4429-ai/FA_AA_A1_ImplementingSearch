#include <divsufsort.h>
#include <filesystem>
#include <vector>
#include <string>
#include <algorithm>
#include <chrono>
#include <cstdint>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>

// dna5 -> char string
static inline void append_dna5_to_string(std::string &dst, auto const &seq)
{
    size_t old = dst.size();
    dst.resize(old + seq.size());
    for (size_t i = 0; i < seq.size(); ++i)
        dst[old + i] = seqan3::to_char(seq[i]);
}

int main(int argc, char const *const *argv)
{
    seqan3::argument_parser parser{"suffixarray_search", argc, argv, seqan3::update_notifications::off};
    parser.info.author = "SeqAn-Team";
    parser.info.version = "1.0.0";

    std::filesystem::path reference_file{};
    parser.add_option(reference_file, '\0', "reference", "Path to the reference FASTA/FASTQ file.");

    std::filesystem::path query_file{};
    parser.add_option(query_file, '\0', "query", "Path to the query FASTA/FASTQ file.");

    size_t number_of_queries{100};
    parser.add_option(number_of_queries, '\0', "query_ct", "Number of queries; if not enough, queries will be duplicated logically.");

    try
    {
        parser.parse();
    }
    catch (seqan3::argument_parser_error const &e)
    {
        seqan3::debug_stream << "Parsing error: " << e.what() << "\n";
        return EXIT_FAILURE;
    }

    // 1) Load reference: concatenate into ONE string (still same as your current behavior)
    seqan3::sequence_file_input reference_stream{reference_file};
    std::string reference;
    reference.reserve(1'000'000); // small hint; will grow

    for (auto &record : reference_stream)
        append_dna5_to_string(reference, record.sequence());

    if (reference.empty())
    {
        seqan3::debug_stream << "Error: reference file contains no sequences.\n";
        return EXIT_FAILURE;
    }

    // 2) Load queries: store ONLY the base queries (do NOT physically duplicate to query_ct)
    seqan3::sequence_file_input query_stream{query_file};
    std::vector<std::string> queries;
    queries.reserve(1024);

    for (auto &record : query_stream)
    {
        std::string q;
        q.reserve(record.sequence().size());
        append_dna5_to_string(q, record.sequence());
        queries.push_back(std::move(q));
    }

    if (queries.empty())
    {
        seqan3::debug_stream << "Error: query file contains no sequences.\n";
        return EXIT_FAILURE;
    }

    // 3) Build suffix array
    seqan3::debug_stream << "Building Suffix Array...\n";
    auto build_start = std::chrono::high_resolution_clock::now();

    std::vector<saidx_t> sa(reference.size());
    auto const *str = reinterpret_cast<sauchar_t const *>(reference.data());
    if (divsufsort(str, sa.data(), static_cast<saidx_t>(reference.size())) != 0)
    {
        seqan3::debug_stream << "Error: divsufsort failed.\n";
        return EXIT_FAILURE;
    }

    auto build_end = std::chrono::high_resolution_clock::now();
    seqan3::debug_stream << "Index Construction time: "
                         << std::chrono::duration<double>(build_end - build_start).count()
                         << " seconds.\n";

    // 4) Search
    size_t total_matches = 0;
    seqan3::debug_stream << "Starting Suffix Array Search for " << number_of_queries << " queries...\n";
    auto search_start = std::chrono::high_resolution_clock::now();

    auto cmp_suffix_vs_pat = [&](saidx_t pos, std::string const &pat)
    {
        size_t m = pat.size();
        size_t n = reference.size();
        size_t end = std::min(n, static_cast<size_t>(pos) + m);
        return std::lexicographical_compare(reference.begin() + pos, reference.begin() + end,
                                            pat.begin(), pat.end());
    };

    auto cmp_pat_vs_suffix = [&](std::string const &pat, saidx_t pos)
    {
        size_t m = pat.size();
        size_t n = reference.size();
        size_t end = std::min(n, static_cast<size_t>(pos) + m);
        return std::lexicographical_compare(pat.begin(), pat.end(),
                                            reference.begin() + pos, reference.begin() + end);
    };

    for (size_t i = 0; i < number_of_queries; ++i)
    {
        std::string const &q = queries[i % queries.size()];
        if (q.empty())
            continue;

        auto it_low = std::lower_bound(sa.begin(), sa.end(), q,
                                       [&](saidx_t sa_pos, std::string const &pat)
                                       {
                                           return cmp_suffix_vs_pat(sa_pos, pat);
                                       });

        auto it_high = std::upper_bound(it_low, sa.end(), q,
                                        [&](std::string const &pat, saidx_t sa_pos)
                                        {
                                            return cmp_pat_vs_suffix(pat, sa_pos);
                                        });

        total_matches += static_cast<size_t>(std::distance(it_low, it_high));
    }

    auto search_end = std::chrono::high_resolution_clock::now();
    seqan3::debug_stream << "Search Finished!\n";
    seqan3::debug_stream << "Total Matches Found: " << total_matches << "\n";
    seqan3::debug_stream << "Search time: "
                         << std::chrono::duration<double>(search_end - search_start).count()
                         << " seconds.\n";

    return 0;
}
