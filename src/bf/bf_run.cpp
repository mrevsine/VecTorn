#include <iostream>
#include <string>
#include <fstream>

#include "../../lib/bloom_filter.hpp"
#include "./common.hpp"

bool use_ref_bf = false;

int main(int argc, char* argv[])
{
	std::ios_base::sync_with_stdio(false);

   	if (argc < 6) {
    	std::cerr << "Usage: " << argv[0] << " <FASTQ> <BF_FILE> <K> <HIT_THRESHOLD> <RESULT_FILE> <OPTIONAL:REF_BF>" << std::endl;
    	return 1;
   	}

	std::string reads_file = argv[1];
   	std::string vec_bf_file = argv[2];
	size_t k = std::stoi(argv[3]);
	size_t vec_threshold = std::stoi(argv[4]);
	std::string result_file = argv[5];
	std::string ref_bf_file;
	if (argc > 6) {
		std::cout << "Using reference bloom filter" << std::endl;
		use_ref_bf = true;
		ref_bf_file = argv[6];
   	}

	bloom_filter vec_bf;
	std::ifstream in_vec_bf(vec_bf_file);
	if (!in_vec_bf.is_open()) {
      	std::cerr << "Error opening " << vec_bf_file << std::endl;
      	return 1;
   	}
	vec_bf.load(in_vec_bf);
	in_vec_bf.close();

	bloom_filter ref_bf;
	if (use_ref_bf) {
		std::ifstream in_ref_bf(ref_bf_file);
		if (!in_ref_bf.is_open()) {
			std::cerr << "Error opening " << ref_bf_file << std::endl;
			return 1;
		}
		ref_bf.load(in_ref_bf);
		in_ref_bf.close();
	}
	
   	std::ifstream in_reads(reads_file);
   	if (!in_reads.is_open()) {
      	std::cerr << "Error opening " << reads_file << std::endl;
      	return 1;
   	}

	std::ofstream out_result(result_file);
   	if (!out_result.is_open()) {
      	std::cerr << "Error opening " << result_file << std::endl;
      	return 1;
   	}

   	std::string read_name;
	std::string read_seq;
	size_t num_lines = 0;
	size_t num_vec_class = 0;
   	while (getline(in_reads, read_name)) {
		// take just the read identifier
		// std::size_t stop_pos = read_name.find(' ');
		// std::string_view read_id = (stop_pos != std::string::npos) ? std::string_view(read_name.data(), stop_pos) : std::string_view(read_name);

		getline(in_reads, read_seq);
		// skip the + and quality lines
		in_reads.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		in_reads.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

		// TODO: threads for kmers
		size_t vec_hit = 0;
		for (size_t i = 0; i < read_seq.length() - k; ++i) {
			std::string kmer = canonical(read_seq.substr(i, k));
			bool hit = vec_bf.contains(kmer);

			if (use_ref_bf) hit = hit && !ref_bf.contains(kmer);
			vec_hit += hit;

			if (vec_hit >= vec_threshold) {
				// out_result << read_id << std::endl;
				out_result << "1" << std::endl;
				++num_vec_class;
				break;
			}
		}
		if (vec_hit < vec_threshold) {
			out_result << "0" << std::endl;
		}
		++num_lines;
   	}
   	in_reads.close();
	out_result.close();

	std::cout << "Finished, vector classified reads found in " << result_file << std::endl;
	std::cout << "Called " << std::to_string(num_vec_class) << "/" << std::to_string(num_lines);
	std::cout << " = " << std::to_string((1.0*num_vec_class/num_lines)*100) << "% Vector Contaminated" << std::endl;

   	return 0;
}