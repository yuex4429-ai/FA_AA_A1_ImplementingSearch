#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <fstream>
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

int main(int argc, char const *const *argv)
{
    seqan3::argument_parser parser{"fmindex_search", argc, argv, seqan3::update_notifications::off};
    parser.info.author = "SeqAn-Team";
    parser.info.version = "1.0.0";

    std::filesystem::path index_path{};
    parser.add_option(index_path, '\0', "index", "Path to the serialized FM-index (.bin).");

    std::filesystem::path query_file{};
    parser.add_option(query_file, '\0', "query", "Path to the query FASTA/FASTQ file.");

    size_t number_of_queries{100};
    parser.add_option(number_of_queries, '\0', "query_ct",
                      "Number of queries; if not enough queries, they will be duplicated.");

    uint8_t number_of_errors{0};
    parser.add_option(number_of_errors, '\0', "errors",
                      "Maximum allowed Hamming errors (mismatches / substitutions only).");

    try
    {
        parser.parse();
    }
    catch (seqan3::argument_parser_error const &ext)
    {
        seqan3::debug_stream << "Parsing error: " << ext.what() << "\n";
        return EXIT_FAILURE;
    }

    // -----------------------------
    // Load queries into unified container
    // -----------------------------
    seqan3::sequence_file_input query_stream{query_file};
    dna5_collection queries{};

    for (auto &record : query_stream)
    {
        auto const &seq = record.sequence();
        queries.emplace_back(seq.begin(), seq.end()); // copy into dna5_vec
    }

    if (queries.empty() && number_of_queries > 0)
    {
        seqan3::debug_stream << "Error: query file contains no sequences.\n";
        return EXIT_FAILURE;
    }

    // Duplicate queries until we have query_ct (if query_ct == 0 -> run 0 searches)
    if (number_of_queries == 0)
    {
        queries.clear();
    }
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

    // -----------------------------
    // Load index (type must match fmindex_construct)
    // -----------------------------
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

    // -----------------------------
    // Search (Hamming-only: substitutions <= k; no insertions/deletions)
    // -----------------------------
    seqan3::configuration const cfg =
        seqan3::search_cfg::max_error_substitution{seqan3::search_cfg::error_count{number_of_errors}} |
        seqan3::search_cfg::max_error_insertion{seqan3::search_cfg::error_count{0}} |
        seqan3::search_cfg::max_error_deletion{seqan3::search_cfg::error_count{0}};

    size_t hits = 0;
    auto t0 = std::chrono::high_resolution_clock::now();

    for (auto const &q : queries)
        for (auto &&res : seqan3::search(q, index, cfg))
        {
            (void)res;
            ++hits;
        }

    auto t1 = std::chrono::high_resolution_clock::now();

    double search_s = std::chrono::duration<double>(t1 - t0).count();
    seqan3::debug_stream << "Search time: " << search_s << " seconds.\n";
    seqan3::debug_stream << "queries=" << queries.size()
                         << " errors=" << static_cast<int>(number_of_errors)
                         << " hits=" << hits << "\n";

    return EXIT_SUCCESS;
}

