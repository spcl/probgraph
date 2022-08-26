#include <iostream>
#include <cassert>
#include <vector>
#include <sys/time.h>
#include "MurmurHash3.h"
#include "sets.hpp"
#include "set_util.hpp"

using namespace std;

// void printSet(StdSet<int> set) {
//   cout << "Printing set:" << endl;
//   for (auto element : set) {
//     cout << element << endl;
//   }
// }

void printHash(int key, int seed) {
  int hash;
  MurmurHash3_x86_32(&key, sizeof(key), seed, &hash);
  cout << "key: " << key << ", hash: " << hash << endl;
}

void printVector(vector<int> &vec) {
  cout << "Printing vector..." << endl << "\t";
  for (auto & e : vec) {
    cout << e << " ";
  }
  cout << endl;
}

template <typename T>
void printSet(KMVSet<T> * set) {
    cout << "Printing KMVSet..." << endl;
    cout << "\tseed: " << set->seed << endl;
    cout << "\telements count: " << set->size() << endl;
    cout << "\tk: " << set->k << endl;
    cout << "\tmin_count: " << set->min_count << endl;
    cout << "\tvalues in min set: ";
    for (int i=0; i<set->min_count; i++) {
        cout << set->values[i] << " ";
    }
    cout << endl;
}

template <typename T>
int exactIntersectionCount(vector<T> &a, vector<T> &b) {
  std::vector<T> i_sect;
  std::set_intersection(a.begin(), a.end(),
                        b.begin(), b.end(),
                        std::back_inserter(i_sect));   
  return i_sect.size();
}

float jaccardEstimate(KMVSet<int> * set_a, KMVSet<int> * set_b, size_t k) {
    std::vector<int> v1(set_a->values, set_a->values + set_a->min_count);
    std::vector<int> v2(set_b->values, set_b->values + set_b->min_count);
    std::vector<int> i_sect;
    std::set_intersection(v1.begin(), v1.end(),
                          v2.begin(), v2.end(),
                          std::back_inserter(i_sect));
    // cout << "intersect_count: " << i_sect.size() << endl << endl;

    float jacc_coeff = (float) i_sect.size() / k;
    return jacc_coeff*(set_a->size() + set_b->size())/(1+jacc_coeff);
}

void testCreation() {
    cout << "Testing KMV creation..." << endl;
    // struct timeval time;
    // gettimeofday(&time, NULL);
    // int seed = time.tv_usec;

    uint32_t seed = 649135;
    // With this seed, the values of MurmurHash3_x86_32 are
    // Hash(1) = 1209617678
    // Hash(2) = -506443777
    // Hash(3) = -289976115
    // Hash(4) = 742609994
    // Hash(5) = 1463702654

    std::vector<int> vec = {1, 2, 3, 4, 5};
    size_t k = 2;

    // The KMinHashSet constructor requires pointers instead of iterators
    int* begin = &vec[0];
    int* end = &vec[0] + vec.size();

    KMVSet<int> * test_set = new KMVSet<int>(begin, end, k, seed);
    assert(test_set->seed == seed);
    assert(test_set->size() == 5);
    assert(test_set->k == k);
    assert(test_set->min_count == (int) k);

    // The k/min_count hash should be in the following form
    assert(test_set->values[0] == -506443777);
    assert(test_set->values[1] == -289976115);
    cout << "Success!" << endl;
    delete test_set;
}

void testSIMD() {
  cout << "Testing intersect_simd4x_count..." << endl;
  std::vector<int> vec_a = {1, 2, 3, 4, 5, 6, 7, 8};
  std::vector<int> vec_b = {3, 6, 8, 10, 13, 56};
  assert(intersect_simd4x_count(&vec_a[0], vec_a.size(), &vec_b[0], vec_b.size()) == 3);

  vec_a = {1, 2, 3, 4, 5, 6, 7, 8};
  vec_b = {3};
  assert(intersect_simd4x_count(&vec_a[0], vec_a.size(), &vec_b[0], vec_b.size()) == 1);

  vec_a = {10};
  vec_b = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  assert(intersect_simd4x_count(&vec_a[0], vec_a.size(), &vec_b[0], vec_b.size()) == 0);
  cout << "Success!" << endl;
}

void testIntersectionCount(vector<int> &vec_a, vector<int> &vec_b, size_t k) {
    cout << "Testing KMV intersection..." << endl;
    int seed = 649135;
    // The KMV constructor requires pointers instead of iterators
    int* begin_a = &vec_a[0];
    int* end_a = &vec_a[0] + vec_a.size();
    KMVSet<int> * set_a = new KMVSet<int>(begin_a, end_a, k, seed);
    int* begin_b = &vec_b[0];
    int* end_b = &vec_b[0] + vec_b.size();
    KMVSet<int> * set_b = new KMVSet<int>(begin_b, end_b, k, seed);
//     printVector(vec_a);
//     printSet(set_a);
//     printVector(vec_b);
//     printSet(set_b);

    int exact_count = exactIntersectionCount(vec_a, vec_b);
    // size_t min = std::min(set_a->min_count, set_b->min_count);
    float jacc_estimate = jaccardEstimate(set_a, set_b, k);
    float actual_estimate = set_a->intersect_count(*set_b);
    cout << "exact count: " << exact_count << endl;
    cout << "desired estimate: " << jacc_estimate << endl;
    cout << "actual estimate: " << actual_estimate << endl;
    assert(jacc_estimate == actual_estimate);
    cout << "Success!" << endl;
    delete set_a;
    delete set_b;
}

void testIntersectionCountMultiple() {
    cout << "Testing KMV multiple intersection..." << endl;
    int seed = 649135;
    size_t k = 14;

    vector<int> vec_a = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14};
    vector<int> vec_b = {3, 4, 5, 6, 7, 8};
    vector<int> vec_c = {1, 2, 3, 4, 5, 6, 7};

    // The KMV constructor requires pointers instead of iterators
    int* begin_a = &vec_a[0];
    int* end_a = &vec_a[0] + vec_a.size();
    KMVSet<int> * set_a = new KMVSet<int>(begin_a, end_a, k, seed);

    int* begin_b = &vec_b[0];
    int* end_b = &vec_b[0] + vec_b.size();
    KMVSet<int> * set_b = new KMVSet<int>(begin_b, end_b, k, seed);

    int* begin_c = &vec_c[0];
    int* end_c = &vec_c[0] + vec_c.size();
    KMVSet<int> * set_c = new KMVSet<int>(begin_c, end_c, k, seed);

    vector<vector<int> *> vecs;
    vecs.push_back(&vec_a);
    vecs.push_back(&vec_b);
    vecs.push_back(&vec_c);
    size_t max_deg = set_util::max_vector_size(vecs);
    int exact_count = set_util::intersect_count_multiple(vecs);
    float expected_estimate = exact_count * max_deg / k;

    // vector<KMVSet<int> *> kminsets = {&set_b, &set_c};
    vector<KMVSet<int> *> kmvsets;
    kmvsets.push_back(set_b);
    kmvsets.push_back(set_c);
    float actual_estimate = set_a->intersect_count_multiple(kmvsets);

    cout << "exact count: " << exact_count << endl;
    cout << "expected estimate: " << expected_estimate << endl;
    cout << "actual estimate: " << actual_estimate << endl;
    // assert(expected_estimate == actual_estimate);
    // cout << "Success!" << endl;
    delete set_a;
    delete set_b;
    delete set_c;
}


int main(){
    testCreation();
    testSIMD();

    std::vector<int> vec_a = {1, 2, 3, 4, 5, 6, 7, 8};
    std::vector<int> vec_b = {3, 6};
    testIntersectionCount(vec_a, vec_b, 3);

    // vec_a = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20};
    // vec_b = {1, 14, 15, 20, 10123, 121231, 12312124, 13345234, 456345, 3456356};
    vec_a = {1, 2, 3, 4, 5};
    vec_b = {1, 10, 12};
    testIntersectionCount(vec_a, vec_b, 5);

    std::vector<int> vec_c;
    vec_c.reserve(100);
    std::vector<int> vec_d;
    vec_d.reserve(100);
    for (int i=0; i<1000; i++) {
      vec_c.push_back(i);
    }
    for (int i=500; i<1500; i++) {
      vec_d.push_back(i);
    }
    testIntersectionCount(vec_c, vec_d, 1000);

    testIntersectionCountMultiple();
    return 0;
}
