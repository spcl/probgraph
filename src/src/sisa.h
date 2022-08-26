#pragma once

#include <assert.h>
#include <algorithm>
#include <utility>
#include <vector>
#include <iterator>
#include <cstdio>
#include <cstring>
#include "benchmark.h"

__always_inline
    NodeID
    UV2E(NodeID u, NodeID v, NodeID n) {
    return u * n + v;
}

__always_inline
    NodeID
    E2U(NodeID E, NodeID n) {
    return E / n;
}

__always_inline
    NodeID
    E2V(NodeID E, NodeID n) {
    return E % n;
}

// A definition of a SISA set (simple model)
class SSSet
{
   public:
    SSSet() = default;
    explicit SSSet(std::vector<NodeID> inputVector): vector(std::move(inputVector))
    {

    }
    SSSet(NodeID* begin, NodeID* end): vector(begin, end)
    {

    }
    explicit SSSet(size_t size): vector(size)
    {

    }
    SSSet(SSSet&& other) noexcept : vector(std::move(other.vector))
    {

    }

    size_t size() const
    {
        return this->vector.size();
    }

    std::vector<NodeID>::iterator begin()
    {
        return this->vector.begin();
    }
    std::vector<NodeID>::const_iterator begin() const
    {
        return this->vector.cbegin();
    }
    std::vector<NodeID>::iterator end()
    {
        return this->vector.end();
    }
    std::vector<NodeID>::const_iterator end() const
    {
        return this->vector.cend();
    }

    SSSet& AddVertex(NodeID v)
    {
        vector.push_back(v);
        return *this;
    }

    SSSet& SubtractVertex(NodeID v)
    {
        vector.erase(std::find(vector.begin(), vector.end(), v));
        return *this;
    }

    void Intersect_Inplace(const SSSet& other)
    {
        //SSSet isect(std::min(set_a.size(), set_b.size()));
        //auto ptr = isect.data();
        auto it = std::set_intersection(this->begin(), this->end(), other.begin(), other.end(), this->begin());
        //assert(ptr == isect.data());
        this->resize(it - this->begin());
       // return isect;
    }

    void append(const SSSet& other)
    {
        vector.insert(vector.end(), other.vector.begin(), other.vector.end());
    }
    void append(NodeID* begin, NodeID* end)
    {
        vector.insert(vector.end(), begin, end);
    }

    template <typename Iterator>
    void erase(Iterator begin, Iterator end)
    {
        vector.erase(begin, end);
    }

    void resize(size_t size)
    {
        vector.resize(size);
    }

    NodeID* data()
    {
        return this->vector.data();
    }

    const NodeID* data() const
    {
        return this->vector.data();
    }

    NodeID& operator[](int i)
    {
        return this->data()[i];
    }

private:
    std::vector<NodeID> vector;
};

template <typename A>
SSSet SISA_Clone(const A& set_a)
{
    SSSet clone(set_a.size());
    std::copy(set_a.begin(), set_a.end(), clone.begin());
    return clone;
}

//STD-based set intersection
template <typename A, typename B>
SSSet SISA_Intersect(const A& set_a, const B& set_b)
{
    SSSet isect(std::min(set_a.size(), set_b.size()));
    auto ptr = isect.data();
    auto it = std::set_intersection(set_a.begin(), set_a.end(), set_b.begin(), set_b.end(), isect.begin());
    assert(ptr == isect.data());
    isect.resize(it - isect.begin());
    return isect;
}

//set intersection size based on the intersection
template <typename A, typename B>
size_t SISA_Intersect_size(const A& set_a, const B& set_b)
{
    return SISA_Intersect(set_a, set_b).size();
}

//STD-based set difference
SSSet SISA_Diff(const SSSet& set_a, const SSSet& set_b)
{
    SSSet diff(std::max(set_a.size(), set_b.size()));
    auto ptr = diff.data();
    auto it = std::set_difference(set_a.begin(), set_a.end(), set_b.begin(), set_b.end(), diff.begin());
    assert(ptr == diff.data());
    diff.resize(it - diff.begin());
    return diff;
}

//STD-based set union
SSSet SISA_Union(const SSSet& set_a, const SSSet& set_b)
{
    SSSet s_union(set_a.size() + set_b.size());
    auto ptr = s_union.data();
    auto it = std::set_union(set_a.begin(), set_a.end(), set_b.begin(), set_b.end(), s_union.begin());
    assert(ptr == s_union.data());
    s_union.resize(it - s_union.begin());
    return s_union;
}

class SISA_Graph {
   public:
    explicit SISA_Graph(const Graph &gapbs_graph)
    {
        this->num_nodes_ = gapbs_graph.num_nodes();

        this->neigh_ids_vector_.reserve(this->num_nodes_);
        this->pred_ids_vector_.reserve(this->num_nodes_);
        this->succ_ids_vector_.reserve(this->num_nodes_);
        this->td.reserve(this->num_nodes_);
        this->n_neighs.reserve(this->num_nodes_);

        for (int64_t u = 0; u < this->num_nodes_; u++)
        {
            SSSet neigh(gapbs_graph.out_neigh(u).begin(), gapbs_graph.out_neigh(u).end());
            SSSet succ, pred;
            if (gapbs_graph.directed())
            {
                neigh.append(gapbs_graph.in_neigh(u).begin(), gapbs_graph.in_neigh(u).end());
            }
            std::sort(neigh.begin(), neigh.end());
            int64_t neigh_sz = 0;
            for (auto v : neigh)
            {
                neigh_sz++;
                if (u > v)
                    pred.AddVertex(v);
                else
                    succ.AddVertex(v);
            }
            n_neighs[u] = neigh_sz;
            if (neigh_sz > 0) {
              td[u].resize(neigh_sz, false);
            }
            this->neigh_ids_vector_.push_back(std::move(neigh));
            this->succ_ids_vector_.push_back(std::move(succ));
            this->pred_ids_vector_.push_back(std::move(pred));
        }
    }

    explicit SISA_Graph(const Graph& g, bool directed)
    {
        this->num_nodes_ = g.num_nodes();

        this->neigh_ids_vector_.reserve(this->num_nodes_);
        this->pred_ids_vector_.reserve(this->num_nodes_);
        this->succ_ids_vector_.reserve(this->num_nodes_);
        
        this->td.reserve(this->num_nodes_);
        this->n_neighs.reserve(this->num_nodes_);

        for(int64_t u = 0; u < this->num_nodes_; u++)
        {
            SSSet neigh(g.out_neigh(u).begin(), g.out_neigh(u).end());
            SSSet outNeigh(g.out_neigh(u).begin(), g.out_neigh(u).end());
            neigh.append( g.in_neigh(u).begin(), g.in_neigh(u).end());
            
            
            
            std::sort(neigh.begin(), neigh.end());
            std::sort(outNeigh.begin(), outNeigh.end());
            SSSet succ, pred;

            int64_t neigh_sz = 0;
            for(auto node : neigh)
            {
                neigh_sz++;
                if(std::binary_search(outNeigh.begin(), outNeigh.end(), node))
                {
                    succ.AddVertex(node);
                }
                else 
                {
                    pred.AddVertex(node);
                }
            }
            n_neighs[u] = neigh_sz;
            if (neigh_sz > 0)
              td[u].resize(neigh_sz, false);

            this->neigh_ids_vector_.push_back(std::move(neigh));
            this->succ_ids_vector_.push_back(std::move(succ));
            this->pred_ids_vector_.push_back(std::move(pred));
            
        }
    }

    NodeID num_nodes() {
        return num_nodes_;
    }

    SSSet& pred(NodeID node)
    {
        return this->pred_ids_vector_[node];
    }

    SSSet& succ(NodeID node)
    {
        return this->succ_ids_vector_[node];
    }

    SSSet& neigh(NodeID node)
    {
        return this->neigh_ids_vector_[node];
    }

    int64_t out_degree(NodeID node)
    {
        return this->pred(node).size();
    }

    SSSet GetEdges()
    {
        NodeID n = this->num_nodes_;
        SSSet edgeSet;
        for (NodeID u = 0; u < n; u++)
        {
            for (auto v : neigh(u))
            {
                edgeSet.AddVertex(UV2E(u, v, n));
            }
        }
        return edgeSet;
    }

    void AddEdge(NodeID u, NodeID v)
    {
        neigh(u).AddVertex(v);
    }

    void RemoveEdge(NodeID u, NodeID v, bool checkIfExists = false)
    {
        //std::cout << neigh(u).size() << std::endl;
        /*
        for (int i = 0; i < neigh(u).size(); i++)
          std::cout << neigh(u)[i] << std::endl;
        */
        neigh(u).erase(std::remove(neigh(u).begin(), neigh(u).end(), v), neigh(u).end());
    }

    void markEdgesforRemoval(NodeID u, NodeID v)
    {
      //assert(v < td[u].size());
      auto v_i = std::lower_bound(neigh(u).begin(), neigh(u).end(), v);
      int64_t vi_ = v_i - neigh(u).begin();
      //std::cout << u << " " << u << " " << v << "vi_: " << vi_ << std::endl;
      td[u][vi_] = true;
    }
    
   std::vector<std::vector<bool> > td; // mark edges to delete
   std::vector<int64_t> n_neighs;

   //private:
    std::vector<SSSet> neigh_ids_vector_, pred_ids_vector_, succ_ids_vector_;
    int64_t num_nodes_;
};

template <typename Iter>
std::string iter2str(Iter begin, Iter end, const std::string &sep = ", ")
{
    std::ostringstream result;
    while (begin != end) {
        result << *begin;
        ++begin;
        if (begin != end)
            result << sep;
    }
    return result.str();
}
std::string set2str(SSSet& set, const std::string &sep = ", ")
{
    return iter2str(set.begin(), set.end(), sep);
}
