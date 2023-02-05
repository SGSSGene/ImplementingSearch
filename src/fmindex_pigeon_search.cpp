#include <sstream>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/search/search.hpp>

int main(int argc, char const* const* argv) {
    seqan3::argument_parser parser{"fmindex_pigeon_search", argc, argv, seqan3::update_notifications::off};

    parser.info.author = "SeqAn-Team";
    parser.info.version = "1.0.0";

    auto index_path = std::filesystem::path{};
    parser.add_option(index_path, '\0', "index", "path to the query file");

    auto reference_file = std::filesystem::path{};
    parser.add_option(reference_file, '\0', "reference", "path to the reference file");

    auto query_file = std::filesystem::path{};
    parser.add_option(query_file, '\0', "query", "path to the query file");

    auto number_of_queries = size_t{100};
    parser.add_option(number_of_queries, '\0', "query_ct", "number of query, if not enough queries, these will be duplicated");

    auto number_of_errors = uint8_t{0};
    parser.add_option(number_of_errors, '\0', "errors", "number of allowed hamming distance errors");

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
    std::vector<std::vector<seqan3::dna5>> reference;
    for (auto& record : reference_stream) {
        reference.push_back(record.sequence());
    }

    // read query into memory
    std::vector<std::vector<seqan3::dna5>> queries;
    for (auto& record : query_stream) {
        queries.push_back(record.sequence());
    }

    // loading fm-index into memory
    using Index = decltype(seqan3::fm_index{std::vector<std::vector<seqan3::dna5>>{}}); // Some hack
    Index index; // construct fm-index
    {
        seqan3::debug_stream << "Loading 2FM-Index ... " << std::flush;
        std::ifstream is{index_path, std::ios::binary};
        cereal::BinaryInputArchive iarchive{is};
        iarchive(index);
        seqan3::debug_stream << "done\n";
    }

    // duplicate input until its large enough
    while (queries.size() < number_of_queries) {
        auto old_count = queries.size();
        queries.resize(2 * old_count);
        std::copy_n(queries.begin(), old_count, queries.begin() + old_count);
    }
    queries.resize(number_of_queries); // will reduce the amount of searches

    seqan3::configuration const cfg = seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{0}};

    // Function that verifies that the hamming distance at a certain position is lower than number_of_errors
    auto verifyDistance = [&](size_t refId, size_t qId, size_t refPos) {
        auto const& r = reference[refId];
        auto const& q = queries[qId];
        if (r.size() < refPos + q.size()) return false; // query is longer than the reference
        size_t errorCount{};
        for (size_t i{0}; i < q.size() and errorCount <= number_of_errors; ++i) {
            if (q[i] != r[refPos+i]) errorCount += 1;
        }
        return errorCount <= number_of_errors;
    };

    //!TODO !ImplementMe use the seqan3::search to find a partial error free hit, verify the rest inside the text
    for (size_t qId{0}; qId < queries.size(); ++qId) {
        auto const& q = queries[qId];
        size_t maxParts = number_of_errors + 1ul;
        for (auto partId{0ul}; partId < maxParts; ++partId) {
            auto startOffset = q.size() * partId / maxParts;
            auto part = std::span{q.begin() + startOffset,
                                  q.begin() + q.size() * (partId+1) / maxParts};
            for (auto result : search(part, index, cfg)) {
                if (result.reference_begin_position() < startOffset) continue;
                if (verifyDistance(result.reference_id(), qId, result.reference_begin_position() - startOffset)) {
                    std::cout << "found query " << qId << " in sequence " << result.reference_id() << " at position " << result.reference_begin_position() - startOffset << "\n";
                }
            }
        }
    }

    return 0;
}
