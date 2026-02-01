#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include <divsufsort.h>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/io/sequence_file/all.hpp>

static inline std::string dna5_to_string(std::vector<seqan3::dna5> const &v)
{
    std::string s;
    s.resize(v.size());
    for (size_t i = 0; i < v.size(); ++i)
        s[i] = seqan3::to_char(v[i]);
    return s;
}

// Write SA as binary: [uint64_t n][uint32_t sa[0..n-1]]
static inline void write_sa(std::filesystem::path const &out, std::vector<uint32_t> const &sa)
{
    std::ofstream os(out, std::ios::binary);
    if (!os)
        throw std::runtime_error("Cannot open output index file.");

    uint64_t n = static_cast<uint64_t>(sa.size());
    os.write(reinterpret_cast<char const *>(&n), sizeof(n));
    os.write(reinterpret_cast<char const *>(sa.data()), sizeof(uint32_t) * sa.size());

    if (!os)
        throw std::runtime_error("Error while writing index file.");
}

int main(int argc, char const *const *argv)
{
    seqan3::argument_parser parser{"suffixarray_construct", argc, argv, seqan3::update_notifications::off};
    parser.info.author = "SeqAn-Team";
    parser.info.version = "1.0.0";

    std::filesystem::path reference_file{};
    parser.add_option(reference_file, '\0', "reference", "Path to reference FASTA/FASTQ (.gz ok)");

    std::filesystem::path index_path{};
    parser.add_option(index_path, '\0', "index", "Path to write suffix array index (.bin)");

    try
    {
        parser.parse();
    }
    catch (seqan3::argument_parser_error const &ext)
    {
        std::cerr << "Parsing error: " << ext.what() << "\n";
        return EXIT_FAILURE;
    }

    // Load reference and concatenate records safely:
    // Insert a separator between contigs to prevent cross-contig false matches.
    // Use '%' as separator (not in DNA alphabet), and '$' as a global sentinel (smallest char).
    seqan3::sequence_file_input reference_stream{reference_file};
    std::string S;
    bool first = true;

    for (auto &rec : reference_stream)
    {
        if (!first)
            S.push_back('%'); // contig separator (queries never contain it)
        first = false;

        S += dna5_to_string(rec.sequence());
    }

    if (S.empty())
    {
        std::cerr << "Error: reference file contains no sequences.\n";
        return EXIT_FAILURE;
    }

    // Global sentinel (assumed smallest, as in lecture/PPT)
    S.push_back('$');

    // Safety: SA stored as uint32_t, require n < 2^32
    if (S.size() > static_cast<size_t>(std::numeric_limits<uint32_t>::max()))
    {
        std::cerr << "Error: reference too long for uint32 suffix array (n >= 2^32).\n";
        return EXIT_FAILURE;
    }

    // Build SA using divsufsort
    std::vector<int32_t> sa32(S.size(), 0);

    auto t0 = std::chrono::high_resolution_clock::now();
    int rc = divsufsort(reinterpret_cast<unsigned char const *>(S.data()), sa32.data(), static_cast<int>(S.size()));
    auto t1 = std::chrono::high_resolution_clock::now();

    if (rc != 0)
    {
        std::cerr << "divsufsort failed with rc=" << rc << "\n";
        return EXIT_FAILURE;
    }

    // Convert to uint32_t for compact storage
    std::vector<uint32_t> sa(S.size());
    for (size_t i = 0; i < sa.size(); ++i)
        sa[i] = static_cast<uint32_t>(sa32[i]);

    try
    {
        write_sa(index_path, sa);
    }
    catch (std::exception const &e)
    {
        std::cerr << "Error writing index: " << e.what() << "\n";
        return EXIT_FAILURE;
    }

    std::chrono::duration<double> dt = t1 - t0;
    std::cout << "Index Construction time: " << dt.count() << " seconds.\n";
    return 0;
}

