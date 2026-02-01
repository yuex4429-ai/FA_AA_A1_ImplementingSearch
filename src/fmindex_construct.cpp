#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>

#include <cereal/archives/binary.hpp>

// -------- unified "dna5_vector" container type (portable) --------
using dna5_vec = std::vector<seqan3::dna5>;
using dna5_collection = std::vector<dna5_vec>;
// ----------------------------------------------------------------

int main(int argc, char const *const *argv)
{
    seqan3::argument_parser parser{"fmindex_construct", argc, argv, seqan3::update_notifications::off};
    parser.info.author = "SeqAn-Team";
    parser.info.version = "1.0.0";

    std::filesystem::path reference_file{};
    parser.add_option(reference_file, '\0', "reference", "Path to the reference FASTA/FASTQ file.");

    std::filesystem::path index_path{};
    parser.add_option(index_path, '\0', "index", "Path to write the serialized FM-index (.bin).");

    try
    {
        parser.parse();
    }
    catch (seqan3::argument_parser_error const &ext)
    {
        seqan3::debug_stream << "Parsing error: " << ext.what() << "\n";
        return EXIT_FAILURE;
    }

    // load reference
    seqan3::sequence_file_input reference_stream{reference_file};

    dna5_collection reference{};
    for (auto &record : reference_stream)
    {
        // record.sequence() could be different container type across versions;
        // copy into our unified dna5_vec.
        auto const &seq = record.sequence();
        reference.emplace_back(seq.begin(), seq.end());
    }

    if (reference.empty())
    {
        seqan3::debug_stream << "Error: reference file contains no sequences.\n";
        return EXIT_FAILURE;
    }

    // build index
    auto t0 = std::chrono::high_resolution_clock::now();
    seqan3::fm_index index{reference};
    auto t1 = std::chrono::high_resolution_clock::now();

    double index_s = std::chrono::duration<double>(t1 - t0).count();
    seqan3::debug_stream << "Index Construction time: " << index_s << " seconds.\n";

    // save index
    try
    {
        std::ofstream os{index_path, std::ios::binary};
        if (!os)
        {
            seqan3::debug_stream << "Error: cannot open index output file: " << index_path << "\n";
            return EXIT_FAILURE;
        }
        cereal::BinaryOutputArchive oarchive{os};
        oarchive(index);
    }
    catch (std::exception const &e)
    {
        seqan3::debug_stream << "Error while saving index: " << e.what() << "\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

