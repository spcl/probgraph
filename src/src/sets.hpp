#include <algorithm>
#include <boost/dynamic_bitset.hpp>
#include <math.h>
#include <vector>
#include "util.hpp"
#include "set_operation.hpp"
#include <functional>
#include <queue>
#include <stddef.h>
#include "MurmurHash3.h"
#include <iterator>
#include <ctime>
#include <set>
#include "set_util.hpp"
#include <cstdlib>

#pragma once

#define LOG_OPS 1 // 0:don't log, 1:do log
namespace set{
//type for object that is used to identify a set that was inserted
template <typename Setid, typename Setelement>
class Sets{
    public:
        Sets(){}
        virtual ~Sets(){}
        virtual Setid addSet(Setelement* begin, Setelement* end) = 0;
//        virtual bool contains(Setid id, Setelement& elem) const  = 0;

//        virtual Setid intersect(Setid first, Setid second) = 0;
        virtual int intersect_count(Setid first, Setid second) const = 0;
        virtual int intersect_count(Setid first, Setid second, Setid third) const = 0;
        virtual Setid intersect(Setid first, Setid second) const =0;

//        virtual Setid merge(Setid first, Setid second) = 0;
//        virtual int merge_count(Setid first, Setid second) = 0;

};

template <typename Setelement>
class Set{
    public:
        virtual int size() const = 0;
        virtual ~Set(){}
        //virtual size_t intersect_count(Set &that) const = 0;
        //virtual Set& intersect(Set&) const = 0;
};
}


template <typename T>
class StdSet: public set::Set<T> {
    public:
        typedef typename std::set<T>::iterator iterator;
        typedef typename std::vector<T>::iterator vector_iterator;
        // This seemingly useless typedef is required by std::set, do not remove
        typedef T value_type;

        StdSet<T>(vector_iterator begin, vector_iterator end) {
            std::set<T> set_from_vector(begin, end);
            this->internal = set_from_vector;
        }

        StdSet<T>(iterator begin, iterator end) {
            std::set<T> set_from_set(begin, end);
            this->internal = set_from_set;
        }

        StdSet<T>(T* begin, T* end) {
            for(T* it = begin; it!=end; ++it){
               this->internal.insert(*it);
            }
            //printf("nv: %.0f m: %d \t nelements: %d estimated: %f\n", nv, m, element_count, nv*(1-pow(1-(float)bitarray.count()/m, (float) m/nv)));

        }

        StdSet<T>() {}

        ~StdSet() {}

        iterator begin() {
            return this->internal.begin();
        }

        iterator end() {
            return this->internal.end();
        }

        int size() const {
            return this->internal.size();
        }

        bool empty() const {
            return this->internal.empty();
        }

        bool contains(T element) const {
            return this->internal.find(element) != this->internal.end();
        } 

        // This is required by stl_iterator, do not remove
        iterator insert(iterator hint, const value_type& value ) {
            return this->internal.insert(hint, value);
        }

        void insert(const T& value) {
            this->internal.insert(value);
        }

        void insert(std::vector<T> v) {
            std::copy(v.begin(), v.end(), std::inserter(this->internal, this->internal.end()));
        }

        void erase(const T& value) {
            this->internal.erase(value);
        }

        iterator erase(iterator it) {
            return this->internal.erase(it);
        }

        StdSet<T> intersect(StdSet<T> &otherSet) {
            StdSet<T> result = StdSet();
            std::set_intersection(
                this->internal.begin(),
                this->internal.end(),
                otherSet.begin(),
                otherSet.end(),
                std::inserter(result, result.begin())
            );
            return result;
        }

        size_t intersect_count(StdSet<T> &otherSet)  {
            StdSet<T> intersection = intersect(otherSet);
            return intersection.size();
        }
        

    private:
        std::set<T> internal;
};


template <typename T>
class VectorSet: public set::Set<T> {
    public:
        typedef typename std::vector<T>::iterator iterator;
        typedef typename std::vector<T>::iterator vector_iterator;
        // This seemingly useless typedef is required by std::set, do not remove
        typedef T value_type;

        VectorSet<T>(T* begin, T* end) {
            for(T* it = begin; it!=end; ++it){
               this->internal.push_back(*it);
            }
        }

        VectorSet<T>() {}

        ~VectorSet() {}

        iterator begin() {
            return this->internal.begin();
        }

        iterator end() {
            return this->internal.end();
        }

        int size() const {
            return this->internal.size();
        }

        bool empty() const {
            return this->internal.empty();
        }

        bool contains(T element) const {
            return this->internal.find(element) != this->internal.end();
        } 

        // This is required by stl_iterator, do not remove
        // iterator insert(iterator hint, const value_type& value ) {
        //     return this->internal.insert(hint, value);
        // }

        iterator insert_sorted(T const& value) {
            return this->internal.insert(std::upper_bound(this->internal.begin(), this->internal.end(), value), value);
        }

        void erase(const T& value) {
            this->internal.erase(value);
        }

        iterator erase(iterator it) {
            return this->internal.erase(it);
        }

        VectorSet<T> intersect(VectorSet<T> &otherSet) {
            VectorSet<T> result = VectorSet();
            std::set_intersection(
                this->internal.begin(),
                this->internal.end(),
                otherSet.begin(),
                otherSet.end(),
                std::inserter(result, result.begin())
            );
            return result;
        }

        size_t intersect_count(VectorSet<T> &otherSet)  {
            VectorSet<T> intersection = intersect(otherSet);
            return intersection.size();
        }

    private:
        std::vector<T> internal;
};

// template <typename se>
// class BPSet: public set::Set<se> {
//     public:
//     class Iterator : public std::iterator<std::input_iterator_tag, se>
//     {
//         public:
//         int base_idx;
//         int state_idx;
//         PackState state;
//         const BPSet* set_;

//         Iterator(const BPSet* set, int bidx, PackState sidx, PackState state):set_(set), base_idx(bidx), state_idx(sidx), state(state){}
//         Iterator(const Iterator& sit) : base_idx(sit.base_idx), state_idx(sit.state_idx), set_(sit.set_), state(sit.state){}
//         Iterator& operator++(){
//             next();
//             return *this;
//         }
//         Iterator operator++(int){
//             Iterator tmp(*this); 
//             operator++();
//             return tmp;
//         }
//         bool operator==(const Iterator& rhs) const {
//             return set_ == rhs.set_ && base_idx == rhs.base_idx && (state_idx == rhs.state_idx || (size_t) base_idx >= set_->block_count);
//         }
//         bool operator!=(const Iterator& rhs) const{
//             return !(*this == rhs);
//         }
//         se operator*(){
//             return (set_->pool_base[base_idx] << PACK_SHIFT) + state_idx;
//         }

//         private:
//         void next(){
//             while(state==0 && base_idx <= set_->block_count) {
//                 if(base_idx >= set_->block_idx){
//                     base_idx = set_->block_count;
//                     return;
//                 }
//                 base_idx ++;
//                 state = set_->pool_state[base_idx];
//                 state_idx = -1;
                
//             }
//             int offset = __builtin_ffs(state);
//             state_idx += offset;
//             state = offset==32? 0: (unsigned int)state >> offset;
//          }
//     };
//     Iterator begin() {
//         if (block_count == 0){
//             return Iterator(this, 0, 0, 0);
//         }else{
//             PackState state = *pool_state;
//             int offset = __builtin_ffs(state);
//             return Iterator(this, 0, offset-1,offset==32? 0:(unsigned int) state>>offset);
//         }
//     };
//     Iterator end() {return Iterator(this, block_count, 0, 0);};

//         BPSet<se>(size_t n){
//             size_t n_blocks = n/sizeof(PackState) + 1;
//             align_malloc((void**)&pool_base, 32, sizeof(int) * (n_blocks+3));
//             align_malloc((void**)&pool_state, 32, sizeof(PackState) * (n_blocks+3));
//             for (se v = 0; v<n; v++){
//                 PackBase v_base = (v >> PACK_SHIFT);
//                 PackState v_bit = ((PackState)1<<(v & PACK_MASK));
//                 if(v == 0){
//                     pool_base[block_idx] = v_base;
//                     pool_state[block_idx] = v_bit;
//                 } else {
//                     if(pool_base[block_idx] == v_base)
//                         pool_state[block_idx] |= v_bit;
//                     else{
//                         pool_base[++block_idx] = v_base;
//                         pool_state[block_idx] = v_bit;
//                     }
//                 }

//             }
//             element_count = n;
//             block_count = block_idx+1;
//         }
         
//         BPSet<se>(se* begin, se* end, float tres=5): treshold(tres){
//             align_malloc((void**)&pool_base, 32, sizeof(int) * (end-begin+3));
//             align_malloc((void**)&pool_state, 32, sizeof(PackState) * (end-begin+3));
//             for (se* v = begin; v<end; v++){
//                 PackBase v_base = (*v >> PACK_SHIFT);
//                 PackState v_bit = ((PackState)1<<(*v & PACK_MASK));
//                 if(v == begin){
//                     pool_base[block_idx] = v_base;
//                     pool_state[block_idx] = v_bit;
//                 } else {
//                     if(pool_base[block_idx] == v_base)
//                         pool_state[block_idx] |= v_bit;
//                     else{
//                         pool_base[++block_idx] = v_base;
//                         pool_state[block_idx] = v_bit;
//                     }
//                 }

//             }
//             element_count = end-begin;
//             block_count = end-begin == 0? 0: block_idx+1;

//         }
//         BPSet<se>(int* pool_base, PackState* pool_state, size_t size){
//             this->pool_base = pool_base;
//             this->pool_state = pool_state;
//             this->block_count = size;
//             this->block_idx = size==0?0:size-1;
//         }
//         ~BPSet<se>(){
//             free(pool_state);
//             free(pool_base);
//         }
//         virtual size_t size() const {return element_count;};
//         size_t n_blocks() const{return block_idx +1;};
//         size_t intersect_count(BPSet<se> &that) const{
//             size_t res=0;
//             const BPSet<se>* larger;
//             const BPSet<se>* smaller;
            
//             if(this->block_idx < that.block_idx){
//                 larger = &that;
//                 smaller = this;
//             }else{
//                 larger = this;
//                 smaller = &that;
//             }
//             if(larger->block_count <= treshold * (smaller->block_count)){
 
// #if SIMD_STATE==2
//                 res += bp_intersect_scalar2x_count(smaller->pool_base, smaller->pool_state, smaller->block_count,
//                     larger->pool_base, larger->pool_state, larger->block_count);
// #elif SIMD_STATE == 4
// #if SIMD_MODE == 0
//                 res += bp_intersect_simd4x_count(smaller->pool_base, smaller->pool_state, smaller->block_count,
//                     larger->pool_base, larger->pool_state , larger->block_count);
// #else
//                 res += bp_intersect_filter_simd4x_count(smaller->pool_base, smaller->pool_state, smaller->block_count,
//                     larger->pool_base, larger->pool_state, larger->block_count);

// #endif
// #else
//                 res += bp_intersect_count(smaller->pool_base, smaller->pool_state, smaller->block_count,
//                     larger->pool_base, larger->pool_state, larger->block_count);

// #endif
//             }
//             else{

//                 res += bp_intersect_galloping_simd4x_count(smaller->pool_base, smaller->pool_state, smaller->block_count,
//                     larger->pool_base, larger->pool_state, larger->block_count);

 
//             }
// #if LOG_OPS==1
//             std::cout << "Intersection count of - Set A size: "<<get_size() << " Set B size: "<< that.get_size() << " Intersection size: " << res << "\n" << std::endl;
// #endif
//            return res; 
//         };
//         virtual BPSet<se>* intersect(BPSet<se>& that) const{
//             const BPSet<se>* smaller;
//             const BPSet<se>* larger;
//             if(this->block_idx < that.block_idx){
//                 larger = &that;
//                 smaller = this;
//             }else{
//                 larger = this;
//                 smaller = &that;
//             }
//             int* c_bases;
//             PackState* c_states;
//             align_malloc((void**)&c_bases, 32, sizeof(int) * (smaller->block_count+3));
//             align_malloc((void**)&c_states, 32, sizeof(PackState)* (smaller->block_count+3));
//             size_t c_size=0;

//             if(larger->block_count > smaller->treshold * (smaller->block_count)){
 
// #if SIMD_STATE==2
//                 c_size = bp_intersect_scalar2x(smaller->pool_base, smaller->pool_state, smaller->block_count,
//                     larger->pool_base, larger->pool_state, larger->block_count, 
//                     c_bases, c_states);
// #elif SIMD_STATE == 4
// #if SIMD_MODE == 0
//                 c_size = bp_intersect_simd4x(smaller->pool_base, smaller->pool_state, smaller->block_count,
//                     larger->pool_base, larger->pool_state , larger->block_count,
//                     c_bases, c_states);
// #else
//                 c_size = bp_intersect_filter_simd4x(smaller->pool_base, smaller->pool_state, smaller->block_count,
//                     larger->pool_base, larger->pool_state, larger->block_count,
//                     c_bases, c_states);

// #endif
// #else
//                 c_size = bp_intersect(smaller->pool_base, smaller->pool_state, smaller->block_count,
//                     larger->pool_base, larger->pool_state, larger->block_count,
//                     c_bases, c_states);

// #endif
//             }else{

//                 c_size = bp_intersect_galloping_simd4x(smaller->pool_base, smaller->pool_state, smaller->block_count,
//                     larger->pool_base, larger->pool_state, larger->block_count, 
//                     c_bases, c_states);

 
//             }
//             BPSet* res =  new BPSet(c_bases, c_states, c_size);
// #if LOG_OPS==1
//             std::cout << "Intersection of - Set A size: "<<get_size();
//             std::cout << " Set B size: "<< that.get_size() << " Intersection size: " << res->get_size() << "\n" << std::endl;
// #endif

//             return res; 
//         } 
//        void lazy_remove(se v){                                                                                                                   
//             int base = v>>5;                                                                                                                      
//             size_t index = 0;
//             while(true){
//                 if (index>=block_count)
//                     return;
//                 if (pool_base[index]==base)
//                     break;
//                 ++index;
//             }
//             PackState k_bit = ((PackState)1<<((v-(base<<5)) & PACK_MASK));                                                                    
// //            while(!compare_and_swap(*(pool_state+index), *(pool_state+index), *(pool_state+index) & ~k_bit));                                    
//        }   
//        void remove_empty_blocks(){
//             size_t to_index = 0;
//             while(to_index < block_count){
//                 if(pool_state[to_index] != 0){
//                     to_index++;
//                 }else{
//                     size_t from_index = to_index+1;
//                     while(from_index < block_count && pool_state[from_index]==0 )
//                         from_index++;
//                     if (from_index>=block_count){
//                         break;
//                     }
//                     pool_state[to_index] = pool_state[from_index];
//                     pool_base[to_index] = pool_base[from_index];
//                     to_index++;
//                 }
//             }
//             block_count = to_index;
//             block_idx = block_count==0? 0:block_count-1;
//         }

// #if LOG_OPS==1
//        int get_size() const{
//            size_t res = 0;
//            for(size_t i = 0; i<block_count; ++i){
//                res += __builtin_popcount(pool_state[i]);
//            }
//            return res;

//        }
// #endif



// //        void lazy_remove(se k){
// //            PackBase k_base = (k >> PACK_SHIFT);
// //            PackState k_bit = ((PackState)1<<(k & PACK_MASK));
// //            size_t idx = find_index(k_base);
// //            if(idx == -1) return;
// //            while(!compare_and_swap(*(pool_state+idx),*(pool_state+idx), *(pool_state+idx) & ~k_bit));
// //        }
//     private:
//         int find_index(size_t s){
//             size_t high = block_idx;
//             size_t low = 0;
//             while(high>low){
//               size_t p = (low + high)/2;
//               if(pool_base[p] == s)
//                   return p;
//               else if(pool_base[p] > s){
//                   high = p-1;
//               }else{
//                   low = p+1;
//               }

//             }
//             if(s == high)
//                 return high;
//             else 
//                 return -1;
//         }
//     public:
//         size_t element_count = 0;  // number of elements in the BPSet object
//         size_t block_idx = 0; // Index of last element in pool_base or pool_state
//         size_t block_count = 0;
//         int* pool_base = NULL; // Array of bases
//         PackState *pool_state = NULL; // Array of states, same length as pool_base
//         float treshold = 5; // if size(larger) > treshold * size(smaller) do galopping intersection


// };

template <typename se>
class OneHashSet: public set::Set<se> {
    public:
        OneHashSet<se>(size_t n, size_t k) {}

        OneHashSet<se>(se* begin, se* end, float k, int seed_): k((size_t)k) {
            // std::cout << "k: " << k << std::endl;
            seed = seed_;
            k = (size_t)k;
            values = (int*) malloc(k*sizeof(int));
            min_count = 0;

            // Populate min-queue with (hash(node), node) pairs
            typedef std::pair<int, int> pi;
            std::priority_queue<pi, std::vector<pi>, std::greater<pi>> pq;
            // std::priority_queue<std::pair<int,int>> pq = std::priority_queue<std::pair<int,int>>();
            for (se* it=begin; it!=end; ++it) {
                int hash;
                MurmurHash3_x86_32(it, sizeof(se), seed, &hash);
                pq.push(std::pair<int,int>(hash, *it));
            }
            // Save the k minimum elements
            for (min_count=0; !pq.empty() && min_count<k; ++min_count) {
                values[min_count] = pq.top().second;
                pq.pop();
            }
            min_count = end-begin > k ? k:end-begin;
            element_count = end-begin;
            std::sort(values, values + min_count);
        }

        ~OneHashSet<se>(){
            free(values);
            values = NULL;
        }

        virtual int size() const { return element_count; }

        float intersect_count_simd(OneHashSet<se> &that) const{
            if (that.k != k) {
                throw std::invalid_argument("intersecting OneHashSets with different k");
            }
            size_t count = intersect_simd4x_count(values, min_count, that.values, that.min_count);
            float j = (float) count / k;
            return j*(element_count + that.element_count)/(1+j);
        };

        float intersect_count(OneHashSet<se> &that) const{
            if (that.k != k) {
                throw std::invalid_argument("intersecting OneHashSets with different k");
            }

            if (min_count == 0 || that.min_count == 0) return 0.0;

            size_t count = 0;
            auto first_this = values;
            auto last_this = values + min_count;
            auto first_that = that.values;
            auto last_that = that.values + that.min_count;
            while (first_this != last_this && first_that != last_that) {
                if (*first_this < *first_that) {
                    ++first_this;
                } else  {
                    if (!(*first_that < *first_this)) {
                        count++;
                    }
                    ++first_that;
                }
            }
            float j = (float) count * 1.0 / std::min(min_count, that.min_count);
            return j*(element_count + that.element_count) * 1.0 /(1+j);
        };

        inline std::vector<int>* get_vector_ptr(int* start, int size) {
            auto res = new std::vector<int>(start, start + size);
            return res;
        }

        float intersect_count_multiple(std::vector<OneHashSet<se> *> &sets, size_t *counter=nullptr) {
            // Convert values in kminhashsets into std::vector
            std::vector<std::vector<int> *> vectorsets;
            vectorsets.reserve(sets.size() + 1); // leave room for this set as well
            for (auto set : sets) {
                auto vec_ptr = get_vector_ptr(set->values, set->min_count);
                vectorsets.push_back(vec_ptr);
            }
            auto this_vec_ptr = get_vector_ptr(this->values, this->min_count);
            vectorsets.push_back(this_vec_ptr);
            size_t max_size = set_util::max_vector_size(vectorsets);
            size_t values_isect_count = set_util::intersect_count_multiple(vectorsets, counter);
            for (auto vec : vectorsets) {
                delete vec;
            }
            return values_isect_count * (float) max_size * 1.0 / (float) k;
        }

        float intersect_count(std::vector<int> &vec, size_t size, size_t *counter=nullptr) {
            auto this_vec_ptr = get_vector_ptr(this->values, this->min_count);
            size_t values_isect_count = set_util::intersect_count(vec, *this_vec_ptr, counter);
            size_t max_size = std::max(size, this_vec_ptr->size());
            delete this_vec_ptr;
            return this->intersect_count_estimate(values_isect_count, max_size);
        }

        // Return the intersection of the values of kminhashsets
        std::pair<std::vector<int>, size_t> intersect_intermediate(std::vector<OneHashSet<se> *> &sets, size_t *counter=nullptr) {
            // Convert values in kminhashsets into std::vector
            std::vector<std::vector<int> *> vectorsets;
            vectorsets.reserve(sets.size() + 1); // leave room for this set as well
            for (auto set : sets) {
                auto vec_ptr = get_vector_ptr(set->values, set->min_count);
                vectorsets.push_back(vec_ptr);
            }
            auto this_vec_ptr = get_vector_ptr(this->values, this->min_count);
            vectorsets.push_back(this_vec_ptr);
            // size_t max_size = set_util::max_vector_size(vectorsets);
            std::pair<std::vector<int>, size_t> result;
            result.first = set_util::intersect_multiple(vectorsets, counter);
            result.second = set_util::max_vector_size(vectorsets);
            for (auto vec : vectorsets) {
                delete vec;
            }
            return result;
        }

        float intersect_count_estimate(size_t isect_count, size_t max_size) {
            return isect_count * (float) max_size / (float) k;
        }

        float union_count(OneHashSet<se> &that) {
            if (that.k != k) {
                throw std::invalid_argument("intersecting OneHashSets with different k");
            }

            size_t count = 0;
            auto first_this = values;
            auto last_this = values + min_count;
            auto first_that = that.values;
            auto last_that = that.values + that.min_count;
            while (first_this != last_this && first_that != last_that) {
                if (*first_this < *first_that) {
                    ++first_this;
                } else  {
                    if (!(*first_that < *first_this)) {
                        count++;
                    }
                    ++first_that;
                }
            }
            float j = (float) count * 1.0 / std::min(min_count, that.min_count);
            return (element_count + that.element_count) * (1 - j * 1.0 / (1 + j));
        };

        size_t total_size() {
            return k*4 + 4*4;
        }

    public:
        size_t element_count = 0;  // number of elements in the KMinHashSet object
        int min_count = 0;
        int* values;
        size_t k; // if size(larger) > treshold * size(smaller) do galopping intersection
        uint32_t seed = time(0);
};

template <typename se>
class KMinHashSet: public set::Set<se> {
public:
    KMinHashSet<se>(size_t n, size_t k) {}

    KMinHashSet<se>(se* begin, se* end, float k_, int seed_): k((size_t)k) {
        // std::cout << "k: " << k << std::endl;
        k = (size_t)k_;
        values = (int*) malloc(k*sizeof(int));
        min_count = 0;

        std::srand(seed_);

        for (int i = 0; i < k; i++) {
            uint32_t seed = std::rand();
            bool first = true;
            int min_hash;
            int64_t min_value = -1;
            for (se *it = begin; it != end; ++it) {
                int hash;
                MurmurHash3_x86_32(it, sizeof(se), seed, &hash);

                if (first || hash < min_hash) {
                    first = false;
                    min_hash = hash;
                    min_value = *it;
                }
            }
            values[i] = min_value;
        }

        min_count = k;
        element_count = end-begin;
        std::sort(values, values+k);
    }

    ~KMinHashSet<se>(){
        free(values);
        values = NULL;
    }

    virtual int size() const { return element_count; }

    float intersect_count_simd(KMinHashSet<se> &that) const{
        if (that.k != k) {
            throw std::invalid_argument("intersecting KMinHashSets with different k");
        }
        size_t count = intersect_simd4x_count(values, min_count, that.values, that.min_count);
        float j = (float) count / k;
        return j*(element_count + that.element_count)/(1+j);
    };

    float intersect_count(KMinHashSet<se> &that) const{
        if (that.k != k) {
            throw std::invalid_argument("intersecting KMinHashSets with different k");
        }
        size_t count = 0;
        auto first_this = values;
        auto last_this = values + min_count;
        auto first_that = that.values;
        auto last_that = that.values + that.min_count;
        while (first_this != last_this && first_that != last_that) {
            if (*first_this < *first_that) {
                ++first_this;
            } else  {
                if (!(*first_that < *first_this)) {
                    count++;
                }
                ++first_that;
            }
        }
        float j = (float) count * 1.0 / k;
        return j*(element_count + that.element_count) * 1.0 / (1 + j);
    };

    float union_count(KMinHashSet<se> &that) const{
        if (that.k != k) {
            throw std::invalid_argument("intersecting KMinHashSets with different k");
        }
        size_t count = 0;
        auto first_this = values;
        auto last_this = values + min_count;
        auto first_that = that.values;
        auto last_that = that.values + that.min_count;
        while (first_this != last_this && first_that != last_that) {
            if (*first_this < *first_that) {
                ++first_this;
            } else  {
                if (!(*first_that < *first_this)) {
                    count++;
                }
                ++first_that;
            }
        }
        float j = (float) count * 1.0 / k;
        return (element_count + that.element_count) * (1 - j * 1.0 / (1 + j));
    }

    inline std::vector<int>* get_vector_ptr(int* start, int size) {
        auto res = new std::vector<int>(start, start + size);
        return res;
    }

    float intersect_count_multiple(std::vector<KMinHashSet<se> *> &sets, size_t *counter=nullptr) {
        // Convert values in kminhashsets into std::vector
        std::vector<std::vector<int> *> vectorsets;
        vectorsets.reserve(sets.size() + 1); // leave room for this set as well
        for (auto set : sets) {
            auto vec_ptr = get_vector_ptr(set->values, set->min_count);
            vectorsets.push_back(vec_ptr);
        }
        auto this_vec_ptr = get_vector_ptr(this->values, this->min_count);
        vectorsets.push_back(this_vec_ptr);
        size_t max_size = set_util::max_vector_size(vectorsets);
        size_t values_isect_count = set_util::intersect_count_multiple(vectorsets, counter);
        for (auto vec : vectorsets) {
            delete vec;
        }
        return values_isect_count * (float) max_size * 1.0 / (float) k;
    }

    float intersect_count(std::vector<int> &vec, size_t size, size_t *counter=nullptr) {
        auto this_vec_ptr = get_vector_ptr(this->values, this->min_count);
        size_t values_isect_count = set_util::intersect_count(vec, *this_vec_ptr, counter);
        size_t max_size = std::max(size, this_vec_ptr->size());
        delete this_vec_ptr;
        return this->intersect_count_estimate(values_isect_count, max_size);
    }

    // Return the intersection of the values of kminhashsets
    std::pair<std::vector<int>, size_t> intersect_intermediate(std::vector<KMinHashSet<se> *> &sets, size_t *counter=nullptr) {
        // Convert values in kminhashsets into std::vector
        std::vector<std::vector<int> *> vectorsets;
        vectorsets.reserve(sets.size() + 1); // leave room for this set as well
        for (auto set : sets) {
            auto vec_ptr = get_vector_ptr(set->values, set->min_count);
            vectorsets.push_back(vec_ptr);
        }
        auto this_vec_ptr = get_vector_ptr(this->values, this->min_count);
        vectorsets.push_back(this_vec_ptr);
        // size_t max_size = set_util::max_vector_size(vectorsets);
        std::pair<std::vector<int>, size_t> result;
        result.first = set_util::intersect_multiple(vectorsets, counter);
        result.second = set_util::max_vector_size(vectorsets);
        for (auto vec : vectorsets) {
            delete vec;
        }
        return result;
    }

    float intersect_count_estimate(size_t isect_count, size_t max_size) {
        return isect_count * (float) max_size / (float) k;
    }

    size_t total_size() {
        return k*4 + 4*4;
    }

public:
    size_t element_count = 0;  // number of elements in the KMinHashSet object
    int min_count = 0;
    int* values;
    size_t k; // if size(larger) > treshold * size(smaller) do galopping intersection
};

template <typename se>
class BlockMinHashSet: public set::Set<se> {
    public:

        BlockMinHashSet<se>(se* begin, se* end, float kfrac, int b_size, int seed_) {
            seed = seed_;
            element_count = end-begin;
            block_size = b_size;
            k = (int) (kfrac*block_size);
            block_bases = (int*) malloc((end-begin)*sizeof(int));
            block_states = (int*) malloc((end-begin)*k*sizeof(int));
            block_sizes = (int*) malloc((end-begin)*sizeof(int));
            std::priority_queue<std::pair<int,int>> pq;
            //skip empty sets
            if(begin == end) return;
            se* it = begin;
            block_idx = 0;
            int bases_idx = 0;
            while(it != end){
                //In one loop iteration either a block will be filled or we skip it.
                //Add the block_idx to the block_bases
                //hash elements while they belong into the current block
                pq = std::priority_queue<std::pair<int,int>>();
                while(pq.empty()){
                    for(;it!=end && *it < block_size*(bases_idx+1); it++){
                        int hash;
                        MurmurHash3_x86_32(it, sizeof(se), seed, &hash);
                        pq.push(std::pair<int,int>(hash, *it));

                    }
                    bases_idx++;
                    if(it==end) break;
                }
                //add mins to the block as long as pq is not empty, update the block idx and the block_sizes accordingly
                block_sizes[block_idx]=0;
                block_bases[block_idx]=bases_idx-1;
                while(!pq.empty() && block_sizes[block_idx] < k){
                    block_states[(block_idx*k) + block_sizes[block_idx]]=pq.top().second;
                    pq.pop();
                    ++block_sizes[block_idx];
                }
                std::sort(block_states+k*block_idx, block_states+k*block_idx + block_sizes[block_idx]);
                block_idx++;
            }
            block_count = block_idx;
    
        }
        BlockMinHashSet<se>(const BlockMinHashSet &old_set){
            element_count = old_set.element_count;
            block_size = old_set.block_size;
            block_count = old_set.block_count;
            block_idx = old_set.block_idx;
            k = old_set.k;
            block_bases = (int*) malloc(sizeof(int)*block_count);
            memcpy(block_bases, old_set.block_bases, block_count);
            block_states = (int*) malloc(sizeof(int)*block_count * k);
            memcpy(block_states, old_set.block_states, block_count*k);
            block_sizes = (int*) malloc(sizeof(int)*block_count);
            memcpy(block_sizes, old_set.block_sizes, block_count);
        }

        ~BlockMinHashSet<se>(){
            free(block_bases);
            block_bases = NULL;
            free(block_states);
            block_states = NULL;
            free(block_sizes);
            block_sizes = NULL;
        }
        virtual int size() const {return element_count;};
        float intersect_count(BlockMinHashSet<se> &that) const{
            if (that.k != k || that.block_size != block_size)
                throw std::invalid_argument("intersecting BlockMinHashSets with different k");
            int i = 0; 
            int l=0;
            float total = 0; 
            while(i<block_count && l<that.block_count){
                if(block_bases[i] == that.block_bases[l]){
                    if(block_sizes[i]<k && that.block_sizes[l]<k){
                        int count = intersect_simd4x_count(block_states+(i*k), block_sizes[i], that.block_states+(l*k), that.block_sizes[l]);
                        total += count;
                    }else{
                    int min = std::min(block_sizes[i], that.block_sizes[l]);
                    int count = intersect_simd4x_count(block_states+(i*k), min, that.block_states+(l*k), min);
                    float j = (float) count / min;
                    total += j*(element_count + that.element_count)/(1+j);
                    }
                    ++i;
                    ++l;
                }
                else if(block_bases[i]<that.block_bases[l]){
                    ++i;
                }else{
                    ++l;
                }
            }
            return total;
        };
        size_t total_size(){
            return block_count * (k*4 + 4 +4 ) + 4*4; // for each block the minhash, the index and size + other variables
        }



    public:
        int element_count = 0;  // number of elements in the BlockMinHashSet object
        int block_size; //number of elements that belong to the same block
        int block_count = 0;
        int block_idx = 0;
        int* block_bases;
        int* block_states;
        int* block_sizes;
        int k = 0; //number of hash values stored
        uint32_t seed = time(0);


};

template <typename se>
class BloomSet: public set::Set<se> {
    public:
        BloomSet<se>(se* begin, se* end, int m_, int k_, std::vector<uint32_t> seeds_, int nv_) {
            element_count = end-begin;
            m = m_;
            //if(element_count!=0) k = m/element_count*0.69314718056; // m/n*ln(2)
            //if(k<1) k = 1;
            k=k_;
            seeds=seeds_;
            nv = (float) nv_;
            bitarray = boost::dynamic_bitset<>(m);
            for(se* it = begin; it!=end; ++it){
               for(size_t j = 0; j<(size_t)k; j++){
                   size_t hash; // hash is always positive
                   MurmurHash3_x86_32(it, sizeof(se), seeds[j], &hash);
                   bitarray.set(hash%m);
               }
            }
            //printf("nv: %.0f m: %d \t nelements: %d estimated: %f\n", nv, m, element_count, nv*(1-pow(1-(float)bitarray.count()/m, (float) m/nv)));
        }

        BloomSet<se>(se* begin, se* end, int m_) {
            element_count = end-begin;
            m = m_;
            if (element_count!=0) {
                k = m/element_count*0.69314718056; // m/n*ln(2)
            }
            if (k<1) {
                k = 1;
            }
            srand(time(0));
            for(int i = 0; i<k; i++){
                seeds.push_back(rand());
            }
            bitarray = boost::dynamic_bitset<>(m);
            for(se* it = begin; it!=end; ++it){
               for(size_t j = 0; j<(size_t)k; j++){
                   size_t hash;
                   MurmurHash3_x86_32(it, 1, seeds[j], &hash);
                   bitarray.set(hash%m);
               }
            }
        }

        BloomSet<se>(const BloomSet &old_set){
            element_count = old_set.element_count;
            k = old_set.k;
            m = old_set.k;
            bitarray = boost::dynamic_bitset<>(old_set.bitarray);
            seeds = std::vector<uint32_t>(old_set.seeds);
        }

        ~BloomSet<se>() {}

        float intersect_count(BloomSet<se> &that) const{
            if (that.k != k || that.m != m)
                throw std::invalid_argument("intersecting BloomSets with different k or m");
            size_t nbits = (bitarray | that.bitarray).count();
            float ret = element_count + that.element_count;
            ret += (float) m * 1.0/k * log(1-(float)nbits * 1.0/m);
            if (ret<0) {
                // this occurs when nbits = m
                return 0;
            }
            return ret;
        };

        float intersect_count_multiple(std::vector<BloomSet<se> *> &sets, size_t *counter=nullptr) {
            // Convert values in kmvSets into std::vector

            auto bitarray_buf = bitarray;
            for (auto set: sets) {
                bitarray_buf &= set->bitarray;
            }
            size_t nbits = bitarray_buf.count();
            float ret = (float) -m * 1.0/k * log(1-(float)nbits * 1.0/m);
            if (ret<0) {
                // this occurs when nbits = m
                return 0;
            }
            return ret;
        };

        float union_count(BloomSet<se> &that) const{
            if (that.k != k || that.m != m)
                throw std::invalid_argument("intersecting BloomSets with different k or m");
            size_t nbits = (bitarray | that.bitarray).count();
            float ret = (float) -m * 1.0/k * log(1-(float)nbits * 1.0/m);
            return ret;
        };

        virtual int size() const { return element_count; }

        size_t total_size() { return bitarray.size()/8 + 4*4 + k*4; }

    public:
        int element_count = 0;  // number of elements in the BlockMinHashSet object
        int k = 0; // number of hash functions
        int m = 0; // length of bit array
        float nv = 0; // number of nodes in graph
        boost::dynamic_bitset<> bitarray;
        std::vector<uint32_t> seeds;
};

// template <typename se>
// class myBloomSet: public set::Set<se> {
//     public:

//         myBloomSet<se>(se* begin, se* end,int m_, int k_, std::vector<uint32_t> seeds_, int nv_) {
//             element_count = end-begin;
//             m = m_;
//             //if(element_count!=0) k = m/element_count*0.69314718056; // m/n*ln(2)
//             //if(k<1) k = 1;
//             k=k_;
//             seeds=seeds_;
//             nv = (float) nv_;
//             nints = m/32 +1;
//             bitarray = new unsigned int[nints](); 
//             for(se* it = begin; it!=end; ++it){
//                for(size_t j = 0; j<(size_t)k; j++){
//                    uint32_t hash;
//                    MurmurHash3_x86_32(it, sizeof(se), seeds[j], &hash);
//                    hash %= m;
//                    bitarray[hash/32] |= 1UL << hash%32;
//                }
//             }
//             //int nbits = 0;
//             //for(size_t i = 0; i<nints; i++) nbits += __builtin_popcount(bitarray[i]);
//             //printf("nv: %.0f m: %d \t nelements: %d estimated: %f\n", nv, m, element_count, nv/k_*(1-pow(1-(float)nbits/m, (float) m/nv)));
//         }
//         ~myBloomSet<se>(){
//             delete[] bitarray;
//         }
//         virtual int size() const {return element_count;};
//         float intersect_count(myBloomSet<se> &that) const{
//             if (that.k != k || that.m != m)
//                 throw std::invalid_argument("intersecting myBloomSets with different k or m");
//             int nbits = 0;
//             for(unsigned int i = 0; i<nints; i++){
//                 nbits += __builtin_popcount(bitarray[i] | that.bitarray[i]);
//             }
            
//             float ret = element_count + that.element_count;
//             //ret -=  nv/k * (1-pow(1-(float)nbits/m, (float)m*k/nv));
//             ret += m/k *log(1-(float)nbits/m);
//             if (ret<0) return 0;
//             return ret;
//         };


//     public:
//         int element_count = 0;  // number of elements in the BlockMinHashSet object
//         int k = 0; //number of hash functions
//         int m = 0; 
//         float nv;
//         unsigned int * bitarray;
//         unsigned int nints;
//         std::vector<uint32_t> seeds;
// };


template <typename se>
class BlockBloomSet: public set::Set<se> {
    public:

        BlockBloomSet<se>(se* begin, se* end, int k_, int m_, int b_size, std::vector<uint32_t> seeds_) {
            element_count = end-begin;
            block_size = b_size;
            k = k_;
            m= m_;

            seeds = seeds_;

            block_bases = (int*) malloc((end-begin)*sizeof(int));
            block_states = std::vector<boost::dynamic_bitset<>>();
            //skip empty sets
            if(begin == end) return;
            block_count = 0;
            int bases_idx = 0;
            block_element_count.push_back(0);
            for(se* it = begin; it != end; it++){
                if(*it>=bases_idx*block_size || block_count==0){
                    // next element is beyond current block
                    while(*it>=(bases_idx+1)*block_size){
                        bases_idx++;
                    }
                    block_bases[block_count] = bases_idx;
                    block_count++;
                    bases_idx++;
                    block_element_count.push_back(0);
                    block_states.push_back(boost::dynamic_bitset<>(m));
                }
                for(size_t j=0; j<(size_t)k;j++){
                    size_t hash;
                    MurmurHash3_x86_32(it, sizeof(se), seeds[j], &hash);
                    block_states.back().set(hash%m);
                    block_element_count.back()++;
                }
                    
            }
        }

        BlockBloomSet<se>(const BlockBloomSet &old_set){
            element_count = old_set.element_count;
            block_size = old_set.block_size;
            block_count = old_set.block_count;
            block_idx = old_set.block_idx;
            k = old_set.k;
            m = old_set.m;
            block_bases = (int*) malloc(sizeof(int)*block_count);
            memcpy(block_bases, old_set.block_bases, block_count);
            block_states = std::vector<boost::dynamic_bitset<>>(old_set.block_states);
        }

        ~BlockBloomSet<se>(){
            free(block_bases);
            block_bases = NULL;
        }
        virtual int size() const {return element_count;};
        float intersect_count(BlockBloomSet<se> &that) const{
            if (that.k != k || that.block_size != block_size || that.m 
                != that.m)
                throw std::invalid_argument("intersecting BlockBloomSets with different k or m");
            int i = 0; 
            int l=0;
            float total = 0; 
            while(i<block_count && l<that.block_count){
                if(block_bases[i] == that.block_bases[l]){
                    size_t nbits = (block_states[i] | that.block_states[l]).count();
                    float ret = block_element_count[i] + that.block_element_count[l];
                    ret += (float) m/k *log(1-(float)nbits/m);
                    if (ret>0) total +=ret;

                    ++i;
                    ++l;
                }
                else if(block_bases[i]<that.block_bases[l]){
                    ++i;
                }else{
                    ++l;
                }
            }
            return total;
        };

        size_t total_size(){
            return block_count * m/8 + block_count*8 +7*4;
        }
        
    public:
        int element_count = 0;  // number of elements in the BlockBloomSet object
        int block_size; //number of elements that belong to the same block
        int block_count = 0;
        int block_idx = 0;
        int* block_bases;
        std::vector<boost::dynamic_bitset<>> block_states;
        std::vector<int> block_element_count;
        int k = 0; //number of hash values stored
        int m = 0;
        std::vector<uint32_t> seeds;
        uint32_t seed = time(0);


};






template <typename  se>

class KMVSet : public set::Set<se>{
public:
    KMVSet<se>(size_t n, size_t k) {}

    KMVSet<se>(se* begin, se* end, float k, int seed_) : k((size_t)k) {
        seed = seed_;
        k = (size_t)k;
        min_count = 0;

        size_t values_size = (end - begin) > k ? k : (end-begin);
        values = (uint32_t*) malloc(values_size * sizeof(uint32_t));

        std::priority_queue<uint32_t,std::vector<uint32_t>,std::greater<uint32_t>> pq;

        for (se* it=begin; it!=end; ++it) {
            uint32_t hash;
            MurmurHash3_x86_32(it, sizeof(se), seed, &hash);
            pq.push(hash);
        }


        // Save the k minimum hashes
        for (min_count=0; !pq.empty() && min_count<k; ++min_count) {
            values[min_count] = pq.top();
            pq.pop();
        }

        min_count = end-begin > k ? k:end-begin;
        element_count = end-begin;
    }

    ~KMVSet<se>(){
        free(values);
        values = NULL;
    }

    virtual int size() const { return element_count; }

    // Find the number of elements in the intersection using Jaccard coef
    float intersect_count_simd(KMVSet<se> &that) const{
        if (that.k != k) {
            throw std::invalid_argument("intersecting KMV sets with different k");
        }
        size_t count = intersect_simd4x_count(values, min_count, that.values, that.min_count);
        float j = (float) count / k;
        return j*(element_count + that.element_count)/(1+j);
    };



    //Calculate the Jaccard coef by itself and use it as approximating the intersection
    float intersect_count(KMVSet<se> &that) const{
        if (that.k != k) {
            throw std::invalid_argument("intersecting KMV sets with different k");
        }

        if (min_count == 0 || that.min_count == 0) {
            return 0.0;
        }

        uint32_t* this_hash = values;
        uint32_t* that_hash = that.values;
        uint32_t current_last_hash = 0;

        auto left_in_this = min_count;
        auto left_in_that = that.min_count;
        auto min_elements = 0;

        for (int i=0; i<k; i++) {
            if (left_in_this == 0 && left_in_that == 0) break;
            min_elements++;
            if (left_in_this == 0 || (*this_hash > *that_hash && left_in_that != 0)) {
                current_last_hash = *that_hash;
                that_hash++;
                left_in_that--;
            } else if (left_in_that == 0 || (*this_hash < *that_hash)) {
                current_last_hash = *this_hash;
                this_hash++;
                left_in_this--;
            } else {
                current_last_hash = *this_hash;
                this_hash++;
                that_hash++;
                left_in_this--;
                left_in_that--;
            }
        }

        double x_k = (min_count == k ? k - 1 : min_count) * 1.0 / ((uint64_t)values[min_count - 1] * 1.0 / ((uint64_t)2 << 32));
        double y_k = (that.min_count == k ? k - 1 : that.min_count) * 1.0 / ((uint64_t)that.values[that.min_count - 1] * 1.0 / ((uint64_t)2 << 32));
        double xuy_k = (min_elements == k ? k - 1 : min_elements) * 1.0 / ((uint64_t)current_last_hash * 1.0 / ((uint64_t)2 << 32));

        return x_k + y_k - xuy_k;
    };


    float union_count(KMVSet<se> &that) const{
        if (that.k != k) {
            throw std::invalid_argument("intersecting KMV sets with different k");
        }

        uint32_t* this_hash = values;
        uint32_t* that_hash = that.values;
        uint32_t current_last_hash = 0;

        auto left_in_this = min_count;
        auto left_in_that = that.min_count;
        auto min_elements = 0;

        for (int i=0; i<k; i++) {
            if (left_in_this == 0 && left_in_that == 0) break;
            min_elements++;
            if (left_in_this == 0 || (*this_hash > *that_hash && left_in_that != 0)) {
                current_last_hash = *that_hash;
                that_hash++;
                left_in_that--;
            } else if (left_in_that == 0 || (*this_hash < *that_hash)) {
                current_last_hash = *this_hash;
                this_hash++;
                left_in_this--;
            } else {
                current_last_hash = *this_hash;
                this_hash++;
                that_hash++;
                left_in_this--;
                left_in_that--;
            }
        }

        double xuy_k = (min_elements == k ? k - 1 : min_elements) * 1.0 / ((uint64_t)current_last_hash * 1.0 / ((uint64_t)2 << 32));
        return xuy_k;
    };



//    The rest of the functions just copied from KMinHash


    std::vector<uint32_t>* get_vector_ptr(uint32_t* start, int size) {
        auto res = new std::vector<uint32_t>(start, start + size);
        return res;
    }

    float intersect_count_multiple(std::vector<KMVSet<se> *> &sets, size_t *counter=nullptr) {
        // Convert values in kmvSets into std::vector
        std::vector<std::vector<uint32_t> *> vectorsets;
        vectorsets.reserve(sets.size() + 1); // leave room for this set as well
        for (auto set : sets) {
            auto vec_ptr = get_vector_ptr(set->values, set->min_count);
            vectorsets.push_back(vec_ptr);
        }
        auto this_vec_ptr = get_vector_ptr(this->values, this->min_count);
        vectorsets.push_back(this_vec_ptr);
        size_t max_size = set_util::max_vector_size(vectorsets);
        size_t values_isect_count = set_util::intersect_count_multiple(vectorsets, counter);
        for (auto vec : vectorsets) {
            delete vec;
        }
        return values_isect_count * (float) max_size * 1.0 / (float) k;
    };

    float intersect_count(std::vector<int> &vec, size_t size, size_t *counter=nullptr) {
        auto this_vec_ptr = get_vector_ptr(this->values, this->min_count);
        size_t values_isect_count = set_util::intersect_count(vec, *this_vec_ptr, counter);
        size_t max_size = std::max(size, this_vec_ptr->size());
        delete this_vec_ptr;
        return this->intersect_count_estimate(values_isect_count, max_size);
    };

    // Return the intersection of the values of kmvSets
    std::pair<std::vector<int>, size_t> intersect_intermediate(std::vector<KMVSet<se> *> &sets, size_t *counter=nullptr) {
        // Convert values in kmvSets into std::vector
        std::vector<std::vector<int> *> vectorsets;
        vectorsets.reserve(sets.size() + 1); // leave room for this set as well
        for (auto set : sets) {
            auto vec_ptr = get_vector_ptr(set->values, set->min_count);
            vectorsets.push_back(vec_ptr);
        }
        auto this_vec_ptr = get_vector_ptr(this->values, this->min_count);
        vectorsets.push_back(this_vec_ptr);
        // size_t max_size = set_util::max_vector_size(vectorsets);
        std::pair<std::vector<int>, size_t> result;
        result.first = set_util::intersect_multiple(vectorsets, counter);
        result.second = set_util::max_vector_size(vectorsets);
        for (auto vec : vectorsets) {
            delete vec;
        }
        return result;
    }

    float intersect_count_estimate(size_t isect_count, size_t max_size) {
        return isect_count * (float) max_size / (float) k;
    }


    float union_count(std::vector<int> &vec, size_t size, size_t *counter=nullptr) {
        return 0;
    };

    size_t total_size() {
        return k*4 + 4*4;
    }





public:
    size_t element_count = 0;  // number of elements in the KMV object
    int min_count = 0;
    uint32_t* values; // the minimum k
    size_t k; // if size(larger) > treshold * size(smaller) do galopping intersection
    uint32_t seed = time(0);
};




template <typename se>
class BlockKMVSet: public set::Set<se> {
public:

    BlockKMVSet<se>(se* begin, se* end, float kfrac, int b_size, int seed_) {
        seed = seed_;
        element_count = end-begin;
        block_size = b_size;
        k = (int) (kfrac*block_size);
        block_bases = (int*) malloc((end-begin)*sizeof(int));
        block_states = (int*) malloc((end-begin)*k*sizeof(int));
        block_sizes = (int*) malloc((end-begin)*sizeof(int));

        std::priority_queue<int,std::vector<int>,std::greater<int>> pq;

        //skip empty sets
        if(begin == end) return;
        se* it = begin;
        block_idx = 0;
        int bases_idx = 0;
        while(it != end){
            //In one loop iteration either a block will be filled or we skip it.
            //Add the block_idx to the block_bases
            //hash elements while they belong into the current block
            while(pq.empty()){
                for(;it!=end && *it < block_size*(bases_idx+1); it++){
                    int hash;
                    MurmurHash3_x86_32(it, sizeof(se), seed, &hash);
                    pq.push(hash);

                }
                bases_idx++;
                if(it==end) break;
            }
            //add mins to the block as long as pq is not empty, update the block idx and the block_sizes accordingly
            block_sizes[block_idx]=0;
            block_bases[block_idx]=bases_idx-1;
            while(!pq.empty() && block_sizes[block_idx] < k){
                block_states[(block_idx*k) + block_sizes[block_idx]]=pq.top();
                pq.pop();
                ++block_sizes[block_idx];
            }
            block_idx++;
        }
        block_count = block_idx;

    }
    BlockKMVSet<se>(const BlockKMVSet &old_set){
        element_count = old_set.element_count;
        block_size = old_set.block_size;
        block_count = old_set.block_count;
        block_idx = old_set.block_idx;
        k = old_set.k;
        block_bases = (int*) malloc(sizeof(int)*block_count);
        memcpy(block_bases, old_set.block_bases, block_count);
        block_states = (int*) malloc(sizeof(int)*block_count * k);
        memcpy(block_states, old_set.block_states, block_count*k);
        block_sizes = (int*) malloc(sizeof(int)*block_count);
        memcpy(block_sizes, old_set.block_sizes, block_count);
    }

    ~BlockKMVSet<se>(){
        free(block_bases);
        block_bases = NULL;
        free(block_states);
        block_states = NULL;
        free(block_sizes);
        block_sizes = NULL;
    }
    virtual int size() const {return element_count;};
    float intersect_count(BlockKMVSet<se> &that) const{
        if (that.k != k || that.block_size != block_size)
            throw std::invalid_argument("intersecting BlockKMVSets with different k");
        int i = 0;
        int l=0;
        float total = 0;
        while(i<block_count && l<that.block_count){
            if(block_bases[i] == that.block_bases[l]){
                if(block_sizes[i]<k && that.block_sizes[l]<k){
                    int count = intersect_simd4x_count(block_states+(i*k), block_sizes[i], that.block_states+(l*k), that.block_sizes[l]);
                    total += count;
                }else{
                    int min = std::min(block_sizes[i], that.block_sizes[l]);
                    int count = intersect_simd4x_count(block_states+(i*k), min, that.block_states+(l*k), min);
                    float j = (float) count / min;
                    total += j*(element_count + that.element_count)/(1+j);
                }
                ++i;
                ++l;
            }
            else if(block_bases[i]<that.block_bases[l]){
                ++i;
            }else{
                ++l;
            }
        }
        return total;
    };
    size_t total_size(){
        return block_count * (k*4 + 4 +4 ) + 4*4; // for each block the kmv, the index and size + other variables
    }



public:
    int element_count = 0;  // number of elements in the BlockKMVSet object
    int block_size; //number of elements that belong to the same block
    int block_count = 0;
    int block_idx = 0;
    int* block_bases;
    int* block_states;
    int* block_sizes;
    int k = 0; //number of hash values stored
    uint32_t seed = time(0);


};
















template<typename T>
void printSet(StdSet<T> &s) {
    for (const T& el : s) { std::cout << el << ", ";}
}

////class BPSets: set::Sets<mSetid, mSetelement>{
//template <typename mSetid, typename mSetelement> 
//class BPSets: public set::Sets<mSetid, mSetelement> {
//   public :
//        BPSets<mSetid, mSetelement>(){
//            align_malloc((void**)&pool_base, 32, sizeof(int) * PACK_NODE_POOL_SIZE);
//            align_malloc((void**)&pool_state, 32, sizeof(PackState) * PACK_NODE_POOL_SIZE);
//        }
//        BPSets<mSetid, mSetelement>(int tres){
//            align_malloc((void**)&pool_base, 32, sizeof(int) * PACK_NODE_POOL_SIZE);
//            align_malloc((void**)&pool_state, 32, sizeof(PackState) * PACK_NODE_POOL_SIZE);
//            treshold = tres;
//        }   
//        
//        virtual ~BPSets<mSetid, mSetelement>(){
//            free(pool_base);
//            free(pool_state);
//                
//        }
//        virtual mSetid addSet(mSetelement* begin, mSetelement* end){
//            setidx.push_back(++cur_packnode_idx);
//            setsizes.push_back(0);
//            
//            for (mSetelement* v = begin; v<end; v++){
//                PackBase v_base = (*v >> PACK_SHIFT);
//                PackState v_bit = ((PackState)1<<(*v & PACK_MASK));
//                if(v == begin){
//                    pool_base[cur_packnode_idx] = v_base;
//                    pool_state[cur_packnode_idx] = v_bit;
//                    setsizes[size]++;
//                } else {
//                    if(pool_base[cur_packnode_idx] == v_base)
//                        pool_state[cur_packnode_idx] |= v_bit;
//                    else{
//                        setsizes[size]++;
//                        pool_base[++cur_packnode_idx] = v_base;
//                        pool_state[cur_packnode_idx] = v_bit;
//                    }
//                }
//
//            }
//            size++;
//            return setsizes.size()-1;
//        }
//        virtual mSetid intersect(mSetid first, mSetid second) const{
//            int* c_states;
//            int* c_bases;
//            int c_max = setsizes[first] > setsizes[second] ? setsizes[second]:setsizes[first];
//            align_malloc((void**)&c_bases, 32, sizeof(int) * c_max);
//            align_malloc((void**)&c_states, 32, sizeof(PackState) * c_max);
//            bp_intersect_galloping_simd4x(pool_base + setidx[first], pool_state + setidx[first], setsizes[first],
//        pool_base + setidx[second], pool_state + setidx[second], setsizes[second],
//        c_bases, c_states);
//            return 0;
//        }
//        virtual int intersect_count(mSetid first, mSetid second) const{
//            assert(first<setsizes.size() && second < setsizes.size());
//            mSetid smaller;
//            mSetid larger;
//            if(setsizes[first] > setsizes[second]){
//                smaller = second;
//                larger = first;
//            }else{
//                smaller = first;
//                larger = second;
//            }
//            int res = 0;
//            if(setsizes[larger] > treshold * setsizes[smaller]){
// 
//#if SIMD_STATE==2
//            res += bp_intersect_scalar2x_count(pool_base + setidx[smaller], pool_state + setidx[smaller], setsizes[smaller],
//                    pool_base + setidx[larger], pool_state + setidx[larger], setsizes[larger]);
//#elif SIMD_STATE == 4
//#if SIMD_MODE == 0
//            res += bp_intersect_simd4x_count(pool_base + setidx[smaller], pool_state + setidx[smaller], setsizes[smaller],
//                    pool_base + setidx[larger], pool_state + setidx[larger], setsizes[larger]);
//#else
//            res += bp_intersect_filter_simd4x_count(pool_base + setidx[smaller], pool_state + setidx[smaller], setsizes[smaller],
//                    pool_base + setidx[larger], pool_state + setidx[larger], setsizes[larger]);
//
//#endif
//#else
//            res += bp_intersect_count(pool_base + setidx[smaller], pool_state + setidx[smaller], setsizes[smaller],
//                    pool_base + setidx[larger], pool_state + setidx[larger], setsizes[larger]);
//
//#endif
//               
//            }
//            else{
//
//            res += bp_intersect_galloping_simd4x_count(pool_base + setidx[smaller], pool_state + setidx[smaller], setsizes[smaller],
//                    pool_base + setidx[larger], pool_state + setidx[larger], setsizes[larger]);
//
// 
//            }
//
//           return res;
//        }
//        virtual int intersect_count(mSetid first, mSetid second, mSetid third) const{
//            assert(first<setsizes.size() && second < setsizes.size());
//            int res = 0;
//            int* tmp_bases = NULL;
//            int tmp_size;
//            PackState* tmp_states = NULL;
//            int s = setsizes[first]>setsizes[second]?setsizes[first]:setsizes[second];
//            align_malloc((void**)&tmp_bases, 32, sizeof(int)*(s+1));
//            align_malloc((void**)&tmp_states, 32, sizeof(PackState)*(s+1));
// #if SIMD_STATE==2
//            tmp_size = bp_intersect_scalar2x(pool_base + setidx[first], pool_state + setidx[first], setsizes[first],
//                    pool_base + setidx[second], pool_state + setidx[second], setsizes[second],
//                    tmp_bases, tmp_states);
//            res += bp_intersect_scalar2x_count(tmp_bases, tmp_states, tmp_size,
//                    pool_base + setidx[third], pool_state + setidx[third], setsizes[third]);
//#elif SIMD_STATE == 4
//#if SIMD_MODE == 0
//            tmp_size = bp_intersect_simd4x(pool_base + setidx[first], pool_state + setidx[first], setsizes[first],
//                    pool_base + setidx[second], pool_state + setidx[second], setsizes[second],
//                    tmp_bases, tmp_states);
//            res += bp_intersect_simd4x_count(tmp_bases, tmp_states, tmp_size,
//                    pool_base + setidx[third], pool_state + setidx[third], setsizes[third]);            
//#else
//            tmp_size = bp_intersect_filter_simd4x(pool_base + setidx[first], pool_state + setidx[first], setsizes[first],
//                    pool_base + setidx[second], pool_state + setidx[second], setsizes[second],
//                    tmp_bases, tmp_states);
//            res += bp_intersect_filter_simd4x_count(tmp_bases, tmp_states, tmp_size,
//                    pool_base + setidx[third], pool_state + setidx[third], setsizes[third]);            
//
//#endif
//#else
//            tmp_size = bp_intersect(pool_base + setidx[first], pool_state + setidx[first], setsizes[first],
//                    pool_base + setidx[second], pool_state + setidx[second], setsizes[second],
//                    tmp_bases, tmp_states);
//            res += bp_intersect_count(tmp_bases, tmp_states, tmp_size,
//                    pool_base + setidx[third], pool_state + setidx[third], setsizes[third]);            
//
//
//#endif           
//            free(tmp_bases);
//            free(tmp_states);
//            return res;
//
//        }
//        int size = 0;  // number of elements in the BPSets object
//        int treshold = 15; //if size(first) > treshold * (second) do scalar intersections
//        std::vector<int> setsizes; // Size of set i
//        std::vector<int> setidx; // Offset of elements in pool_base that belong to set i
//        int cur_packnode_idx = -1; // current length of pool_base-1
//        int* pool_base = NULL; // Array of bases
//        PackState *pool_state = NULL; // Array of states, same length as pool_base
//};



