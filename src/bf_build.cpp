#include <iostream>
#include <string>
#include <fstream>

#include "../lib/bloom_filter.hpp"
#include "./common.hpp"

int main(int argc, char* argv[])
{
   if (argc < 5) {
      std::cerr << "Usage: " << argv[0] << " <kmers_file> <out_file> <number_of_elements> <false_pos_probabilitity>" << std::endl;
      return 1;
   }

   std::string kmer_file = argv[1];
   std::string ds_file = argv[2];

	bloom_parameters params;
	params.projected_element_count = std::stoi(argv[3]);
	params.false_positive_probability = std::stod(argv[4]);
   
	params.random_seed = SEED;
	params.compute_optimal_parameters();

	bloom_filter filter(params);
	
   std::ios_base::sync_with_stdio(false);
   std::ifstream in_file(kmer_file);
   if (!in_file.is_open()) {
      std::cerr << "Error opening " <<  std::endl;
      return 1;
   }

   std::string kmer;
   while (getline(in_file, kmer)) {
      filter.insert(canonical(kmer));   
   }
   in_file.close();

   std::cout << "Filter size (bytes): " << filter.serialize(NullStream::ns) << std::endl;

   std::ofstream out_file(ds_file);
   filter.serialize(out_file);
   out_file.close();

   return 0;
}