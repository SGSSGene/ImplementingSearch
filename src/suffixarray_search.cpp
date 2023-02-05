#include <divsufsort.h>
#include <sstream>
#include <span>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/search/search.hpp>

int main(int argc, char const* const* argv) {
    seqan3::argument_parser parser{"suffixarray_search", argc, argv, seqan3::update_notifications::off};

    parser.info.author = "SeqAn-Team";
    parser.info.version = "1.0.0";

    auto reference_file = std::filesystem::path{};
    parser.add_option(reference_file, '\0', "reference", "path to the reference file");

    auto query_file = std::filesystem::path{};
    parser.add_option(query_file, '\0', "query", "path to the query file");

    auto number_of_queries = size_t{100};
    parser.add_option(number_of_queries, '\0', "query_ct", "number of query, if not enough queries, these will be duplicated");

    try {
         parser.parse();
    } catch (seqan3::argument_parser_error const& ext) {
        seqan3::debug_stream << "Parsing error. " << ext.what() << "\n";
        return EXIT_FAILURE;
    }

    // loading our files
    auto reference_stream = seqan3::sequence_file_input{reference_file};
    auto query_stream     = seqan3::sequence_file_input{query_file};

    // read reference into memory
    // Attention: we are concatenating all sequences into one big combined sequence
    //            this is done to simplify the implementation of suffix_arrays
    std::vector<seqan3::dna5> reference;
    for (auto& record : reference_stream) {
        auto r = record.sequence();
        reference.insert(reference.end(), r.begin(), r.end());
    }

    // read query into memory
    std::vector<std::vector<seqan3::dna5>> queries;
    for (auto& record : query_stream) {
        queries.push_back(record.sequence());
    }

    // duplicate input until its large enough
    while (queries.size() < number_of_queries) {
        auto old_count = queries.size();
        queries.resize(2 * old_count);
        std::copy_n(queries.begin(), old_count, queries.begin() + old_count);
    }
    queries.resize(number_of_queries); // will reduce the amount of searches

    // Array that should hold the future suffix array
    std::vector<saidx_t> suffixarray;
    suffixarray.resize(reference.size()); // resizing the array, so it can hold the complete SA

    //!TODO !ImplementMe implement suffix array sort
    //Hint, if can use libdivsufsort (already integrated in this repo)
    //      https://github.com/y-256/libdivsufsort
    //      To make the `reference` compatible with libdivsufsort you can simply
    //      cast it by calling:
    sauchar_t const* str = reinterpret_cast<sauchar_t const*>(reference.data());
    divsufsort(str, suffixarray.data(), suffixarray.size());


    for (size_t qId{0}; qId < queries.size(); ++qId) {
        // Only compares the first n-th elements
        auto comparisonFunction = [&](std::span<seqan3::dna5 const> lhs, std::span<seqan3::dna5 const> rhs) {
            for (size_t i{0}; i < std::min(lhs.size(), rhs.size()); ++i) {
                if (lhs[i] > rhs[i]) return false;
                if (lhs[i] < rhs[i]) return true;
            }
            return false;
        };
        // Transforms the suffixarray integer, into a string view
        auto projectionFunction = [&](saidx_t pos) -> std::span<seqan3::dna5 const> {
            return {reference.begin() + pos, reference.end()};
        };

        auto rangeOnSA = std::ranges::equal_range(suffixarray, std::span{queries[qId]}, comparisonFunction, projectionFunction);

        for (auto pos : rangeOnSA) {
            std::cout << "found query " << qId << " in sequence " << 0 << " at position " << pos << "\n";
        }
    }

    return 0;
}
