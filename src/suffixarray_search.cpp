#include <sstream>
#include <divsufsort.h>

#include <seqan3/std/filesystem>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/search/search.hpp>

struct my_dna4 : seqan3::sequence_file_input_default_traits_dna {
    using sequence_alphabet = seqan3::dna4;
};

// prints out all occurences of query inside of ref
void findOccurences(std::vector<seqan3::dna4> const& ref, std::vector<seqan3::dna4> const& query) {
    std::string ref_str;
    ref_str.reserve(ref.size()); // This is not required, but speeds things up
    for (auto c : ref) {         // Loop over the vector and convert each entry to a char
        ref_str += c.to_char();
    }
//    for (auto 
    //!TODO Implement me
    for (int i{0}; i < ref.size() - query.size(); ++i) {
        bool isQueryFound{true};
        for (int j{0}; j < query.size(); ++j) {
            if (ref[i+j] != query[j]) {
                isQueryFound = false;
                break;
            }
        }
        if (isQueryFound) {
            std::cout << "found at : " << i << "\n";
        }
    }
}

int main(int, char**) {

    // paths to our files
    auto reference_file = std::filesystem::path{"../data/text.dna4.short.fasta.gz"};
    auto query_file     = std::filesystem::path{"../data/sampled_illumina_reads.fasta.gz"};

    // loading our files
    auto reference_stream = seqan3::sequence_file_input<my_dna4>{reference_file};
    auto query_stream     = seqan3::sequence_file_input<my_dna4>{query_file};

    // read reference into memory
    // Attention: we are concatenating all sequences into one big combined sequence
    //            this is done to simplify the implementation of suffix_arrays
    std::vector<seqan3::dna4> reference;
    for (auto& record : reference_stream) {
        auto r = record.sequence();
        reference.insert(reference.end(), r.begin(), r.end());
    }

    // read query into memory
    std::vector<std::vector<seqan3::dna4>> queries;
    for (auto& record : query_stream) {
        queries.push_back(record.sequence());
    }

    //!TODO here adjust the number of searches
    queries.resize(100); // will reduce the amount of searches

    // Construct your Suffix Array
//    std::vector

    //! search for all occurences of queries inside of reference
/*    for (auto& r : reference) {
        for (auto& q : queries) {
            findOccurences(r, q);
        }
    }*/

    return 0;
}
