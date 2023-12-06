#ifndef INCLUDE_COMMON_HPP
#define INCLUDE_COMMON_HPP

#include <iostream>
#include <string>

const unsigned int SEED = 5381;

namespace NullStream {
   // to get file sizes, without writing
   class NullBuffer : public std::streambuf
   {
      public:
      int overflow (int c) { return c; }
   };
   
   NullBuffer nb;
   std::ostream ns(&nb);
}

// DBJHash borrowed from bloom_filter library: www.partow.net/programming/hashfunctions/
unsigned int DJBHash(const char* str, unsigned int length)
{
   unsigned int hash = SEED;
   unsigned int i    = 0;

   for (i = 0; i < length; ++str, ++i)
   {
      hash = ((hash << 5) + hash) + (*str);
   }

   return hash;
}

char complement(char c)
{   
    switch(c)
    {   
      case 'A':
         return 'T';
      case 'T':
         return 'A';
      case 'G':
         return 'C';
      case 'C':
         return 'G';
      default:
         return 'A';
    }   
}

char alpha_upper(char c) {
   return c^0x20;    // xor the 32nd bit, akin to c - 'A' = c - 32
}

// Check that kmer is upper case and canonical
std::string canonical(std::string str) {
	size_t k = str.length();
	std::string upper_case(str);
	std::string rev_comp;
	rev_comp.resize(k);

	for(size_t i = 0; i < str.length(); ++i) {
		char upper = alpha_upper(str[i]);
		upper_case[i] = upper;
		rev_comp[k-i-1] = complement(upper);
	}
	return (DJBHash(upper_case.c_str(), k) <= DJBHash(rev_comp.c_str(), k)) ? upper_case : rev_comp;
}

// Check that kmer is upper case and canonical (string view version)
std::string canonical(std::string_view str) {
	size_t k = str.length();
	std::string upper_case(str);
	std::string rev_comp;
	rev_comp.resize(k);

	for(size_t i = 0; i < str.length(); ++i) {
		char upper = alpha_upper(str[i]);
		upper_case[i] = upper;
		rev_comp[k-i-1] = complement(upper);
	}
	return (DJBHash(upper_case.c_str(), k) <= DJBHash(rev_comp.c_str(), k)) ? upper_case : rev_comp;
}

#endif