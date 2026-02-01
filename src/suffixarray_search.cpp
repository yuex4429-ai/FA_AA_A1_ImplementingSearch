#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <chrono>
#include <cstdint>
#include <algorithm>

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

// Read SA binary: [uint64_t n][uint32_t sa...]
static inline std::vector<uint32_t> read_sa(std::filesystem::path const &in)
{
    std::ifstream is(in, std::ios::binary);
    if (!is)
        throw std::runtime_error("Cannot open index file.");

    uint64_t n = 0;
    is.read(reinterpret_cast<char *>(&n), sizeof(n));
    if (n == 0)
        throw std::runtime_error("Index file corrupt (n=0).");

    std::vector<uint32_t> sa(n);
    is.read(reinterpret_cast<char *>(sa.data()), sizeof(uint32_t) * n);
    return sa;
}

// compare suffix S[pos..] with pattern P (prefix compare)
// -1 suffix < P, 0 match (P is prefix), +1 suffix > P
static inline int cmp_suffix_pattern(std::string const &S, size_t pos, std::string const &P)
{
    size_t n = S.size(), m = P.size();
    size_t i = 0;
    while (i < m && pos + i < n)
    {
        char a = S[pos + i], b = P[i];
        if (a < b)
            return -1;
        if (a > b)
            return +1;
        ++i;
    }
    if (i == m)
        return 0; // P fully matched
    return -1;    // suffix ended
}

// returns [LP, RP] inclusive; empty => RP = -1
static inline std::pair<size_t, size_t> find_interval(std::string const &S,
                                                      std::vector<uint32_t> const &sa,
                                                      std::string const &P)
{
    size_t n = sa.size();

    // LP = first suffix not < P
    size_t L = 0, R = n;
    while (L < R)
    {
        size_t M = (L + R) / 2;
        int c = cmp_suffix_pattern(S, sa[M], P);
        if (c == -1)
            L = M + 1;
        else
            R = M;
    }
    size_t LP = L;

    // first_gt = first suffix > P (treat match as <=)
    L = 0;
    R = n;
    while (L < R)
    {
        size_t M = (L + R) / 2;
        int c = cmp_suffix_pattern(S, sa[M], P);
        if (c == +1)
            R = M;
        else
            L = M + 1;
    }
    size_t first_gt = L;

    if (LP >= first_gt)
        return {0, static_cast<size_t>(-1)};
    return {LP, first_gt - 1};
}

int main(int argc, char const *const *argv)
{
    seqan3::argument_parser parser{"suffixarray_search", argc, argv, seqan3::update_notifications::off};
    parser.info.author = "SeqAn-Team";
    parser.info.version = "1.0.0";

    std::filesystem::path reference_file{};
    parser.add_option(reference_file, '\0', "reference", "Path to reference FASTA/FASTQ (.gz ok)");

    std::filesystem::path index_path{};
    parser.add_option(index_path, '\0', "index", "Path to suffix array index (.bin)");

    std::filesystem::path query_file{};
    parser.add_option(query_file, '\0', "query", "Path to query FASTA/FASTQ (.gz ok)");

    size_t number_of_queries{100};
    parser.add_option(number_of_queries, '\0', "query_ct", "Number of queries; duplicate if needed");

    try
    {
        parser.parse();
    }
    catch (seqan3::argument_parser_error const &ext)
    {
        std::cerr << "Parsing error: " << ext.what() << "\n";
        return EXIT_FAILURE;
    }

    // Load reference (concatenate)
    seqan3::sequence_file_input reference_stream{reference_file};
    std::string S;
    for (auto &rec : reference_stream)
        S += dna5_to_string(rec.sequence());
    S.push_back('$'); // must match construct :contentReference[oaicite:4]{index=4}

    // Load SA
    auto sa = read_sa(index_path);

    // Load queries
    seqan3::sequence_file_input query_stream{query_file};
    std::vector<std::string> queries;
    for (auto &rec : query_stream)
        queries.push_back(dna5_to_string(rec.sequence()));

    // duplicate until enough
    while (queries.size() < number_of_queries)
    {
        size_t old = queries.size();
        queries.resize(2 * old);
        std::copy_n(queries.begin(), old, queries.begin() + old);
    }
    queries.resize(number_of_queries);

    uint64_t total_hits = 0;

    auto t0 = std::chrono::high_resolution_clock::now();
    for (auto const &P : queries)
    {
        auto [LP, RP] = find_interval(S, sa, P);
        if (RP != static_cast<size_t>(-1))
            total_hits += (RP - LP + 1);
    }
    auto t1 = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> dt = t1 - t0;
    std::cout << "Search time: " << dt.count() << " seconds.\n";
    std::cout << "queries=" << queries.size() << " hits=" << total_hits << "\n";
    return 0;
}

