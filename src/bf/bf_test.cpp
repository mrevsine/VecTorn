#include <iostream>
#include <string>
#include <fstream>

#include "../../lib/bloom_filter.hpp"
#include "./common.hpp"

int main(int argc, char* argv[])
{
	if (argc < 2) {
    	std::cerr << "Usage: " << argv[0] << " <kmers_file> <out_file>" << std::endl;
    	return 1;
	}

	std::string kmer_file = argv[1];
	std::string ds_filename = argv[2];
	
	std::ios_base::sync_with_stdio(false);
	std::ifstream ds_file(ds_filename);
	if (!ds_file.is_open()) {
      std::cerr << "Error opening " <<  std::endl;
      return 1;
   	}

	bloom_filter filter;
	filter.load(ds_file);
	ds_file.close();

	std::ifstream in_file(kmer_file);
   	if (!in_file.is_open()) {
      	std::cerr << "Error opening " <<  std::endl;
      	return 1;
   	}

   	std::string kmer;
   	size_t total = 0;
   	size_t correct = 0;
   	while (getline(in_file, kmer)) {
      	correct += filter.contains(canonical(kmer));
	  	++total; 
   	}
   	in_file.close();

   	std::cout << correct << "/" << total << std::endl;
	filter.printValues();

   	return 0;
}