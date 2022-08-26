#include <iostream>
#include <cassert>
#include <vector>
#include <sys/time.h>
#include "MurmurHash3.h"
#include "sets.hpp"


using namespace std;

uint32_t seed_1 = 649135;
// With this seed, the (size_t) cast values of MurmurHash3_x86_32 are
// key: 1,  hash: 140730108038414, %m=4: 2
// key: 2,  hash: 140732686944255, %m=4: 3
// key: 3,  hash: 140732903411917, %m=4: 1
// key: 4,  hash: 140729641030730, %m=4: 2
// key: 5,  hash: 140730362123390, %m=4: 2
// key: 6,  hash: 140731016125830, %m=4: 2
// key: 7,  hash: 140729963672348, %m=4: 0
// key: 8,  hash: 140729301550232, %m=4: 0
// key: 9,  hash: 140729181848190, %m=4: 2
// key: 10, hash: 140732788631270, %m=4: 2
  uint32_t seed_2 = 281384;
// With this seed, the (size_t) cast values of MurmurHash3_x86_32 are
// key: 1,  hash: 140730821808818, %m=4: 2
// key: 2,  hash: 140732521918113, %m=4: 1
// key: 3,  hash: 140730067728711, %m=4: 3
// key: 4,  hash: 140732184101426, %m=4: 2
// key: 5,  hash: 140731753419814, %m=4: 2
// key: 6,  hash: 140731821477518, %m=4: 2
// key: 7,  hash: 140732054749834, %m=4: 2
// key: 8,  hash: 140732290319225, %m=4: 1
// key: 9,  hash: 140733171837782, %m=4: 2
// key: 10, hash: 140731334607166, %m=4: 2


// void printSet(StdSet<int> set) {
//   cout << "Printing set:" << endl;
//   for (auto element : set) {
//     cout << element << endl;
//   }
// }

void printPositiveHash(int key, int seed, int divisor) {
  size_t hash;
  MurmurHash3_x86_32(&key, sizeof(key), seed, &hash);
  cout << "key: " << key << ", hash: " << hash
       << ", %m=" << divisor << ": " << hash%divisor << endl;
}

void printVector(vector<int> &vec) {
  cout << "Printing vector..." << endl << "\t";
  for (auto & e : vec) {
    cout << e << " ";
  }
  cout << endl;
}

void printBitset(boost::dynamic_bitset<> &bitset) {
  cout << "Printing bitset..." << endl << "\t";
  cout << bitset << endl;
}


// template <typename T>
// void printSet(KMinHashSet<T> * set) {
//   cout << "Printing KMinHashSet..." << endl;
//   cout << "\tseed: " << set->seed << endl;
//   cout << "\telements count: " << set->size() << endl;
//   cout << "\tk: " << set->k << endl;
//   cout << "\tmin_count: " << set->min_count << endl;
//   cout << "\tvalues in min set: ";
//   for (int i=0; i<set->min_count; i++) {
//     cout << set->values[i] << " ";
//   }
//   cout << endl;
// }

template <typename T>
int exactIntersectionCount(vector<T> &a, vector<T> &b) {
  std::vector<T> i_sect;
  std::set_intersection(a.begin(), a.end(),
                        b.begin(), b.end(),
                        std::back_inserter(i_sect));   
  return i_sect.size();
}

void testCreation() {
  cout << "Testing BloomSet creation..." << endl;
  std::vector<int> vec = {1, 2, 3, 4, 5};
  int k = 2;
  int m = 4;
  int nv = 4;

  // struct timeval seed;
  // gettimeofday(&seed, NULL);
  // srand(seed.tv_usec);
  // std::vector<uint32_t> seeds;
  // for (int i = 0; i<k; i++) {
  //   seeds.push_back(rand());
  // }
  
  std::vector<uint32_t> seeds;
  seeds.push_back(seed_1);
  seeds.push_back(seed_2);

  // bitvector under seed_1:  1 1 1 0
  // bitvector under seed_2:  1 1 1 0
  // combined bitvector:      1 1 1 0
  boost::dynamic_bitset<> target_bitset(m);
  target_bitset.set(1);
  target_bitset.set(2);
  target_bitset.set(3);

  // The BloomSet constructor requires pointers instead of iterators
  int* begin = &vec[0];
  int* end = &vec[0] + vec.size();
  BloomSet<int> * test_set = new BloomSet<int>(begin, end, m, k, seeds, nv);
  // printBitset(test_set->bitarray);

  assert(test_set->bitarray == target_bitset);
  assert(test_set->seeds[0] == seed_1);
  assert(test_set->seeds[1] == seed_2);
  assert(test_set->size() == 5);
  assert(test_set->k == k);
  assert(test_set->m == m);
  assert(test_set->nv == nv);

  cout << "Success!" << endl;
}

void testIntersectionCount() {
  int k = 2;
  int m = 4;
  int nv = 4;
  std::vector<uint32_t> seeds;
  seeds.push_back(seed_1);
  seeds.push_back(seed_2);

  std::vector<int> vec_1 = {1, 8}; // bitset: 0 1 1 1
  int* begin_1 = &vec_1[0];
  int* end_1 = &vec_1[0] + vec_1.size();
  BloomSet<int> * set_1 = new BloomSet<int>(begin_1, end_1, m, k, seeds, nv);

  // printPositiveHash(1, seed_1, m);
  // printPositiveHash(1, seed_2, m);
  // printPositiveHash(2, seed_1, m);
  // printPositiveHash(2, seed_2, m);
  // printPositiveHash(8, seed_1, m);
  // printPositiveHash(8, seed_2, m);

  std::vector<int> vec_2 = {1, 2};  // bitset: 1 1 1 0
  int* begin_2 = &vec_2[0];
  int* end_2 = &vec_2[0] + vec_2.size();
  BloomSet<int> * set_2 = new BloomSet<int>(begin_2, end_2, m, k, seeds, nv);

  int exact_count = exactIntersectionCount(vec_1, vec_2);

  printBitset(set_1->bitarray);
  printBitset(set_2->bitarray);

  size_t n_bits = (set_1->bitarray | set_2->bitarray).count();
  assert (n_bits == 4);
  cout << "num bits intersection bitsets: " << n_bits << endl;

  float target_estimate = set_1->element_count + set_2->element_count;
  target_estimate += (float) m/k * log(1-(float)n_bits/m);
  if (target_estimate<0) { target_estimate = 0; } // this occurs when nbits = m
  
  float actual_estimate = set_1->intersect_count(*set_2);
  cout << "exact count: " << exact_count << endl;
  cout << "target estimate: " << target_estimate << endl;
  cout << "actual estimate: " << actual_estimate << endl;
  assert(target_estimate == actual_estimate);
  cout << "Success!" << endl;
}

int main(){
    testCreation();
    testIntersectionCount();
    return 0;
}
