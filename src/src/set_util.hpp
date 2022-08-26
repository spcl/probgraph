#include <vector>
#pragma once

#ifndef COUNT
  #define COUNT 0
#endif

namespace set_util {

// Returns the intersection count of two _sorted_ containers
template <typename T>
size_t intersect_count(T &first, T &last, size_t *counter=nullptr) {
  size_t count = 0;
  auto first_this = first.begin();
  auto last_this = first.end();
  auto first_that = last.begin();
  auto last_that = last.end();
  while (first_this != last_this && first_that != last_that) {
    #if COUNT
      if (counter != NULL) {
        *counter = *counter + 1;
      }
    #endif
    if (*first_this < *first_that) {
      ++first_this;
    } else  {
      #if COUNT
        if (counter != NULL) {
          *counter = *counter + 1;
        }
      #endif
      if (!(*first_that < *first_this)) {
        count++;
      }
      ++first_that;
    }
  }
  return count;
}

// Adapted from https://en.cppreference.com/w/cpp/algorithm/set_intersection
template<class InputIt1, class InputIt2, class OutputIt>
void intersect(InputIt1 first1, InputIt1 last1,
               InputIt2 first2, InputIt2 last2,
               OutputIt d_first, size_t *counter=nullptr) {
  while (first1 != last1 && first2 != last2) {
    #if COUNT
      if (counter != NULL) {
        *counter = *counter + 1;
      }
    #endif
    if (*first1 < *first2) {
      ++first1;
    } else  {
      #if COUNT
        if (counter != NULL) {
          *counter = *counter + 1;
        }
      #endif
      if (!(*first2 < *first1)) {
        *d_first++ = *first1++;
      }
      ++first2;
    }
  }
}

typedef std::vector<int> vectorset;

// Return the max size of a list of vectors
inline size_t max_vector_size(std::vector<vectorset *> &sets) {
  size_t max_degree = 0;
  for (auto set_p : sets) {
    max_degree = std::max(set_p->size(), max_degree);
  }
  return max_degree;
}

inline size_t max_vector_size(std::vector<std::vector<uint32_t>*> &sets) {
    size_t max_degree = 0;
    for (auto set_p : sets) {
        max_degree = std::max(set_p->size(), max_degree);
    }
    return max_degree;
}

void print_vvset(std::vector<vectorset *> &sets) {
  std::cout << "printing vector of vectorsets" << std::endl;
  for (auto set : sets) {
    std::cout << "...printing set:" << std::endl;
    for (auto e : (*set)) {
      std::cout << e << std::endl;
    }
  }
  std::cout << "finished printing" << std::endl;
}

// Returns the intersection count of multiple _sorted_ containers
// This could be further space optimized
size_t intersect_count_multiple(std::vector<vectorset *> &sets, size_t *counter=nullptr) {
  size_t max_size = max_vector_size(sets);
  // auto last_isect = *sets[0];
  std::vector<int> last_isect;
  last_isect.reserve(max_size);
  std::copy(sets[0]->begin(), sets[0]->end(), std::back_inserter(last_isect));
  std::vector<int> curr_isect;
  curr_isect.reserve(max_size);

  for (size_t i=1; i<sets.size(); i++) {
    intersect(last_isect.begin(), last_isect.end(),
              sets[i]->begin(), sets[i]->end(),
              std::back_inserter(curr_isect), counter);

    std::swap(last_isect, curr_isect);
    curr_isect.clear();
  }
  return last_isect.size();
}

size_t intersect_count_multiple(std::vector<std::vector<uint32_t>*> &sets, size_t *counter=nullptr) {
    size_t max_size = max_vector_size(sets);
    // auto last_isect = *sets[0];
    std::vector<int> last_isect;
    last_isect.reserve(max_size);
    std::copy(sets[0]->begin(), sets[0]->end(), std::back_inserter(last_isect));
    std::vector<int> curr_isect;
    curr_isect.reserve(max_size);

    for (size_t i=1; i<sets.size(); i++) {
        intersect(last_isect.begin(), last_isect.end(),
                  sets[i]->begin(), sets[i]->end(),
                  std::back_inserter(curr_isect), counter);

        std::swap(last_isect, curr_isect);
        curr_isect.clear();
    }
    return last_isect.size();
}

// Returns the intersection of multiple _sorted_ containers
std::vector<int> intersect_multiple(std::vector<vectorset *> &sets, size_t *counter=nullptr) {
  size_t max_size = max_vector_size(sets);
  // auto last_isect = *sets[0];
  std::vector<int> last_isect;
  last_isect.reserve(max_size);
  std::copy(sets[0]->begin(), sets[0]->end(), std::back_inserter(last_isect));
  std::vector<int> curr_isect;
  curr_isect.reserve(max_size);

  for (size_t i=1; i<sets.size(); i++) {
    intersect(last_isect.begin(), last_isect.end(),
              sets[i]->begin(), sets[i]->end(),
              std::back_inserter(curr_isect), counter);

    std::swap(last_isect, curr_isect);
    curr_isect.clear();
  }
  return last_isect;
}

}
