#include <iostream>
#include <string>
#include <fstream>

#include "../lib/bloom_filter.hpp"
#include "./common.hpp"

bool use_ref_bf = false;

int main(int argc, char* argv[])
{
	std::ios_base::sync_with_stdio(false);

   	if (argc < 5) {
    	std::cerr << "Usage: " << argv[0] << " <FASTQ> <BF_FILE> <K> <HIT_THRESHOLD> <RESULT_FILE> <OPTIONAL:REF_BF>" << std::endl;
    	return 1;
   	}

	std::string reads_file = argv[1];
   	std::string vec_bf_file = argv[2];
	size_t k = std::stoi(argv[3]);
	size_t vec_threshold = std::stoi(argv[4]);
	std::string result_file = argv[5];
	std::string ref_bf_file;
	if (argc <= 6) {
		use_ref_bf;
		ref_bf_file = argv[5];
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
   	while (getline(in_reads, read_name)) {
		// take just the read identifier
		std::size_t stop_pos = read_name.find(' ');
		std::string_view read_id = (stop_pos != std::string::npos) ? std::string_view(read_name.data(), stop_pos) : std::string_view(read_name);

		getline(in_reads, read_seq);
		// skip the + and quality lines
		in_reads.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		in_reads.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

		// TODO: threads for kmers
		size_t vec_hit = 0;
		for (size_t i = 0; i < read_seq.length() - k; ++i) {
			std::string kmer = canonical(std::string_view(read_seq.data() + i, k));
			bool hit = vec_bf.contains(kmer);

			if (use_ref_bf) hit = hit && !ref_bf.contains(kmer);
			vec_hit += hit;

			if (vec_hit >= vec_threshold) {
				out_result << read_id << std::endl;
				break;
			}
		}
   	}
   	in_reads.close();
	out_result.close();

   	return 0;
}