/*********************************************************************
* Software License Agreement (BSD License)
*
*  Copyright (c) 2011, Rice University
*  All rights reserved.
*
*  Redistribution and use in source and binary forms, with or without
*  modification, are permitted provided that the following conditions
*  are met:
*
*   * Redistributions of source code must retain the above copyright
*     notice, this list of conditions and the following disclaimer.
*   * Redistributions in binary form must reproduce the above
*     copyright notice, this list of conditions and the following
*     disclaimer in the documentation and/or other materials provided
*     with the distribution.
*   * Neither the name of the Rice University nor the names of its
*     contributors may be used to endorse or promote products derived
*     from this software without specific prior written permission.
*
*  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
*  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
*  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
*  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
*  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
*  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
*  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
*  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
*  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
*  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
*  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
*  POSSIBILITY OF SUCH DAMAGE.
*********************************************************************/

/* Author: Mark Moll, Bryant Gipson, Shi Shenglei */

#ifndef OMPL_DATASTRUCTURES_NEAREST_NEIGHBORS_GNAT_
#define OMPL_DATASTRUCTURES_NEAREST_NEIGHBORS_GNAT_

#include "ompl/datastructures/NearestNeighbors.h"
#include "ompl/datastructures/GreedyKCenters.h"
#ifdef GNAT_SAMPLER
#include "ompl/datastructures/PDF.h"
#endif
#include <algorithm>
#include <iostream>
#include <queue>
#include <random>
#include <unordered_set>
#include <utility>
#include "ompl/util/Exception.h"

namespace ompl
{
    /** \brief Geometric Near-neighbor Access Tree (GNAT), a data
        structure for nearest neighbor search.

        If GNAT_SAMPLER is defined, then extra code will be enabled to sample
        elements from the GNAT with probability inversely proportial to their
        local density.

        @par External documentation
        S. Brin, Near neighbor search in large metric spaces, in <em>Proc. 21st
        Conf. on Very Large Databases (VLDB)</em>, pp. 574–584, 1995.

        B. Gipson, M. Moll, and L.E. Kavraki, Resolution independent density
         estimation for motion planning in high-dimensional spaces, in
        <em>IEEE Intl. Conf. on Robotics and Automation</em>, 2013.
        [[PDF]](http://www.kavrakilab.org/publications/gipson-moll2013resolution-independent-density.pdf)
    */
    template <typename _T>
    class NearestNeighborsGNAT : public NearestNeighbors<_T>
    {
    protected:
        /// \cond IGNORE
        // internally, we use a priority queue for nearest neighbors, paired
        // with their distance to the query point
        using NearQueue = std::priority_queue<std::pair<double, const _T *>>;

        // another internal data structure is a priority queue of nodes to
        // check next for possible nearest neighbors
        class Node;
        using NodeDist = std::pair<Node *, double>;
        struct NodeDistCompare
        {
            bool operator()(const NodeDist &n0, const NodeDist &n1) const
            {
                return (n0.second - n0.first->maxRadius_) > (n1.second - n1.first->maxRadius_);
            }
        };
        using NodeQueue = std::priority_queue<NodeDist, std::vector<NodeDist>, NodeDistCompare>;
        /// \endcond

    public:
        NearestNeighborsGNAT(unsigned int degree = 8, unsigned int minDegree = 4, unsigned int maxDegree = 12,
                             unsigned int maxNumPtsPerLeaf = 50, unsigned int removedCacheSize = 500,
                             bool rebalancing = false
#ifdef GNAT_SAMPLER
                             ,
                             double estimatedDimension = 6.0
#endif
                             )
          : NearestNeighbors<_T>()
          , degree_(degree)
          , minDegree_(std::min(degree, minDegree))
          , maxDegree_(std::max(maxDegree, degree))
          , maxNumPtsPerLeaf_(maxNumPtsPerLeaf)
          , rebuildSize_(rebalancing ? maxNumPtsPerLeaf * degree : std::numeric_limits<std::size_t>::max())
          , removedCacheSize_(removedCacheSize)
#ifdef GNAT_SAMPLER
          , estimatedDimension_(estimatedDimension)
#endif
        {
        }

        ~NearestNeighborsGNAT() override
        {
            delete tree_;
        }
        /// \brief Set the distance function to use
        void setDistanceFunction(const typename NearestNeighbors<_T>::DistanceFunction &distFun) override
        {
            NearestNeighbors<_T>::setDistanceFunction(distFun);
            pivotSelector_.setDistanceFunction(distFun);
            if (tree_)
                rebuildDataStructure();
        }

        void clear() override
        {
            clearWithoutDisabled();
            disabled_.clear();
        }

        bool reportsSortedResults() const override
        {
            return true;
        }

        void add(const _T &data) override
        {
            if (tree_)
            {
                if (isRemoved(data))
                    rebuildDataStructure();
                tree_->add(*this, data);
            }
            else
            {
                tree_ = new Node(degree_, maxNumPtsPerLeaf_, data);
                size_ = 1;
            }
        }
        void add(const std::vector<_T> &data) override
        {
            if (tree_)
                NearestNeighbors<_T>::add(data);
            else if (!data.empty())
            {
                tree_ = new Node(degree_, maxNumPtsPerLeaf_, data[0]);
#ifdef GNAT_SAMPLER
                tree_->subtreeSize_ = data.size();
#endif
                tree_->data_.insert(tree_->data_.end(), data.begin() + 1, data.end());
                size_ += data.size();
                if (tree_->needToSplit(*this))
                    tree_->split(*this);
            }
        }
        /// \brief Rebuild the internal data structure.
        void rebuildDataStructure()
        {
            for(auto it = disabled_.begin(); it != disabled_.end();)
            {
                if (isRemoved(*(*it)))
                    it = disabled_.erase(it);
                else 
                    it++;
            }
            std::vector<_T> lst;
            list(lst);
            clear();
            add(lst);
        }

        void enable(const _T &data) override
        {
            auto pos = disabled_.find(&data);
            if (pos != disabled_.end())
                disabled_.erase(pos);
        }

        void disable(const _T &data) override
        {
            disabled_.insert(&data);
        }

        /// \brief Remove data from the tree.
        /// The element won't actually be removed immediately, but just marked
        /// for removal in the removed_ cache. When the cache is full, the tree
        /// will be rebuilt and the elements marked for removal will actually
        /// be removed.
        bool remove(const _T &data) override
        {
            if (size_ == 0u)
                return false;
            NearQueue nbhQueue;
            // find data in tree
            bool isPivot = nearestKInternal(data, 1, nbhQueue);
            const _T *d = nbhQueue.top().second;
            if (*d != data)
                return false;
            removed_.insert(d);
            size_--;
            // if we removed a pivot or if the capacity of removed elements
            // has been reached, we rebuild the entire GNAT
            if (isPivot || removed_.size() >= removedCacheSize_)
                rebuildDataStructure();
            return true;
        }

        _T nearest(const _T &data) const override
        {
            if (size_ - disabled_.size())
            {
                NearQueue nbhQueue;
                nearestKInternal(data, 1, nbhQueue);
                if (!nbhQueue.empty())
                    return *nbhQueue.top().second;
            }
            throw Exception("No elements found in nearest neighbors data structure");
        }

        /// Return the k nearest neighbors in sorted order
        void nearestK(const _T &data, std::size_t k, std::vector<_T> &nbh) const override
        {
            nbh.clear();
            if (k == 0)
                return;
            if (size_ - disabled_.size())
            {
                NearQueue nbhQueue;
                nearestKInternal(data, k, nbhQueue);
                postprocessNearest(nbhQueue, nbh);
            }
        }

        /// Return the nearest neighbors within distance \c radius in sorted order
        void nearestR(const _T &data, double radius, std::vector<_T> &nbh) const override
        {
            nbh.clear();
            if (size_ - disabled_.size())
            {
                NearQueue nbhQueue;
                nearestRInternal(data, radius, nbhQueue);
                postprocessNearest(nbhQueue, nbh);
            }
        }

        _T nearestDisabled(const _T &data) const override
        {
            if (!disabled_.empty())
            {
                NearQueue nbhQueue;
                nearestKInternal(data, 1, nbhQueue, true);
                if (!nbhQueue.empty())
                    return *nbhQueue.top().second;
            }
            throw Exception("No elements found in disabled nearest neighbors data structure");
        }

        /// Return the k nearest neighbors in sorted order
        void nearestKDisabled(const _T &data, std::size_t k, std::vector<_T> &nbh) const override
        {
            nbh.clear();
            if (k == 0)
                return;
            if (!disabled_.empty())
            {
                NearQueue nbhQueue;
                nearestKInternal(data, k, nbhQueue, true);
                postprocessNearest(nbhQueue, nbh);
            }
        }

        /// Return the nearest neighbors within distance \c radius in sorted order
        void nearestRDisabled(const _T &data, double radius, std::vector<_T> &nbh) const override
        {
            nbh.clear();
            if (!disabled_.empty())
            {
                NearQueue nbhQueue;
                nearestRInternal(data, radius, nbhQueue, true);
                postprocessNearest(nbhQueue, nbh);
            }
        }

        std::size_t size() const override
        {
            return size_;
        }

        virtual std::size_t disabledSize() const
        {
            return disabled_.size();
        }

#ifdef GNAT_SAMPLER
        /// Sample an element from the GNAT.
        const _T &sample(RNG &rng) const
        {
            if (!size())
                throw Exception("Cannot sample from an empty tree");
            else
                return tree_->sample(*this, rng);
        }
#endif

        void list(std::vector<_T> &data) const override
        {
            data.clear();
            data.reserve(size());
            if (tree_)
                tree_->list(*this, data);
        }

        /// \brief Print a GNAT structure (mostly useful for debugging purposes).
        friend std::ostream &operator<<(std::ostream &out, const NearestNeighborsGNAT<_T> &gnat)
        {
            if (gnat.tree_)
            {
                out << *gnat.tree_;
                if (!gnat.removed_.empty())
                {
                    out << "Elements marked for removal:\n";
                    for (const auto &elt : gnat.removed_)
                        out << *elt << '\t';
                    out << std::endl;
                }
            }
            return out;
        }

        // for debugging purposes
        void integrityCheck()
        {
            std::vector<_T> lst;
            std::unordered_set<const _T *> tmp;
            // get all elements, including those marked for removal
            removed_.swap(tmp);
            list(lst);
            // check if every element marked for removal is also in the tree
            for (const auto &elt : tmp)
            {
                unsigned int i;
                for (i = 0; i < lst.size(); ++i)
                    if (lst[i] == *elt)
                        break;
                if (i == lst.size())
                {
                    // an element marked for removal is not actually in the tree
                    std::cout << "***** FAIL!! ******\n" << *this << '\n';
                    for (const auto &l : lst)
                        std::cout << l << '\t';
                    std::cout << std::endl;
                }
                assert(i != lst.size());
            }
            // restore
            removed_.swap(tmp);
            // get elements in the tree with elements marked for removal purged from the list
            list(lst);
            if (lst.size() != size_)
                std::cout << "#########################################\n" << *this << std::endl;
            assert(lst.size() == size_);
        }

    protected:
        using GNAT = NearestNeighborsGNAT<_T>;

        /// Return true iff data has been marked for removal.
        bool isRemoved(const _T &data) const
        {
            return !removed_.empty() && removed_.find(&data) != removed_.end();
        }

        /// Return true iff data is disabled.
        bool isDisabled(const _T &data) const
        {
            return !disabled_.empty() && disabled_.find(&data) != disabled_.end();
        }

        void clearWithoutDisabled()
        {
            if (tree_)
            {
                delete tree_;
                tree_ = nullptr;
            }
            size_ = 0;
            removed_.clear();
            if (rebuildSize_ != std::numeric_limits<std::size_t>::max())
                rebuildSize_ = maxNumPtsPerLeaf_ * degree_;
        }

        /// \brief Return in nbhQueue the k nearest neighbors of data.
        /// For k=1, return true if the nearest neighbor is a pivot.
        /// (which is important during removal; removing pivots is a
        /// special case).
        bool nearestKInternal(const _T &data, std::size_t k, NearQueue &nbhQueue, bool disabled = false) const
        {
            bool isPivot;
            double dist;
            NodeDist nodeDist;
            NodeQueue nodeQueue;

            dist = NearestNeighbors<_T>::distFun_(data, tree_->pivot_);
            if (disabled == isDisabled(tree_->pivot_))
                isPivot = tree_->insertNeighborK(nbhQueue, k, tree_->pivot_, data, dist);
            tree_->nearestK(*this, data, k, nbhQueue, nodeQueue, isPivot, disabled);
            while (!nodeQueue.empty())
            {
                dist = nbhQueue.top().first;  // note the difference with nearestRInternal
                nodeDist = nodeQueue.top();
                nodeQueue.pop();
                if (nbhQueue.size() == k && (nodeDist.second > nodeDist.first->maxRadius_ + dist ||
                                             nodeDist.second < nodeDist.first->minRadius_ - dist))
                    continue;
                nodeDist.first->nearestK(*this, data, k, nbhQueue, nodeQueue, isPivot, disabled);
            }
            return isPivot;
        }
        /// \brief Return in nbhQueue the elements that are within distance radius of data.
        void nearestRInternal(const _T &data, double radius, NearQueue &nbhQueue, bool disabled = false) const
        {
            double dist = radius;  // note the difference with nearestKInternal
            NodeQueue nodeQueue;
            NodeDist nodeDist;

            if (disabled == isDisabled(tree_->pivot_))
                tree_->insertNeighborR(nbhQueue, radius, tree_->pivot_,
                        NearestNeighbors<_T>::distFun_(data, tree_->pivot_));
            tree_->nearestR(*this, data, radius, nbhQueue, nodeQueue, disabled);
            while (!nodeQueue.empty())
            {
                nodeDist = nodeQueue.top();
                nodeQueue.pop();
                if (nodeDist.second > nodeDist.first->maxRadius_ + dist ||
                    nodeDist.second < nodeDist.first->minRadius_ - dist)
                    continue;
                nodeDist.first->nearestR(*this, data, radius, nbhQueue, nodeQueue, disabled);
            }
        }
        /// \brief Convert the internal data structure used for storing neighbors
        /// to the vector that NearestNeighbor API requires.
        void postprocessNearest(NearQueue &nbhQueue, std::vector<_T> &nbh) const
        {
            nbh.resize(nbhQueue.size());
            for (auto it = nbh.rbegin(); it != nbh.rend(); it++, nbhQueue.pop())
                *it = *nbhQueue.top().second;
        }

        /// The class used internally to define the GNAT.
        class Node
        {
        public:
            /// \brief Construct a node of given degree with at most
            /// \e capacity data elements and with given pivot.
            Node(int degree, int capacity, _T pivot)
              : degree_(degree)
              , pivot_(std::move(pivot))
              , minRadius_(std::numeric_limits<double>::infinity())
              , maxRadius_(-minRadius_)
              , minRange_(degree, minRadius_)
              , maxRange_(degree, maxRadius_)
#ifdef GNAT_SAMPLER
              , subtreeSize_(1)
              , activity_(0)
#endif
            {
                // The "+1" is needed because we add an element before we check whether to split
                data_.reserve(capacity + 1);
            }

            ~Node()
            {
                for (auto &child : children_)
                    delete child;
            }

            /// \brief Update minRadius_ and maxRadius_, given that an element
            /// was added with distance dist to the pivot.
            void updateRadius(double dist)
            {
                if (minRadius_ > dist)
                    minRadius_ = dist;
#ifndef GNAT_SAMPLER
                if (maxRadius_ < dist)
                    maxRadius_ = dist;
#else
                if (maxRadius_ < dist)
                {
                    maxRadius_ = dist;
                    activity_ = 0;
                }
                else
                    activity_ = std::max(-32, activity_ - 1);
#endif
            }
            /// \brief Update minRange_[i] and maxRange_[i], given that an
            /// element was added to the i-th child of the parent that has
            /// distance dist to this Node's pivot.
            void updateRange(unsigned int i, double dist)
            {
                if (minRange_[i] > dist)
                    minRange_[i] = dist;
                if (maxRange_[i] < dist)
                    maxRange_[i] = dist;
            }
            /// Add an element to the tree rooted at this node.
            void add(GNAT &gnat, const _T &data)
            {
#ifdef GNAT_SAMPLER
                subtreeSize_++;
#endif
                if (children_.empty())
                {
                    data_.push_back(data);
                    gnat.size_++;
                    if (needToSplit(gnat))
                    {
                        if (!gnat.removed_.empty())
                            gnat.rebuildDataStructure();
                        else if (gnat.size_ >= gnat.rebuildSize_)
                        {
                            gnat.rebuildSize_ <<= 1;
                            gnat.rebuildDataStructure();
                        }
                        else
                            split(gnat);
                    }
                }
                else
                {
                    std::vector<double> dist(children_.size());
                    double minDist = dist[0] = gnat.distFun_(data, children_[0]->pivot_);
                    int minInd = 0;

                    for (unsigned int i = 1; i < children_.size(); ++i)
                        if ((dist[i] = gnat.distFun_(data, children_[i]->pivot_)) < minDist)
                        {
                            minDist = dist[i];
                            minInd = i;
                        }
                    for (unsigned int i = 0; i < children_.size(); ++i)
                        children_[i]->updateRange(minInd, dist[i]);
                    children_[minInd]->updateRadius(minDist);
                    children_[minInd]->add(gnat, data);
                }
            }
            /// Return true iff the node needs to be split into child nodes.
            bool needToSplit(const GNAT &gnat) const
            {
                unsigned int sz = data_.size();
                return sz > gnat.maxNumPtsPerLeaf_ && sz > degree_;
            }
            /// \brief The split operation finds pivot elements for the child
            /// nodes and moves each data element of this node to the appropriate
            /// child node.
            void split(GNAT &gnat)
            {
                typename GreedyKCenters<_T>::Matrix dists(data_.size(), degree_);
                std::vector<unsigned int> pivots;

                children_.reserve(degree_);
                gnat.pivotSelector_.kcenters(data_, degree_, pivots, dists);
                for (unsigned int &pivot : pivots)
                    children_.push_back(new Node(degree_, gnat.maxNumPtsPerLeaf_, data_[pivot]));
                degree_ = pivots.size();  // in case fewer than degree_ pivots were found
                for (unsigned int j = 0; j < data_.size(); ++j)
                {
                    unsigned int k = 0;
                    for (unsigned int i = 1; i < degree_; ++i)
                        if (dists(j, i) < dists(j, k))
                            k = i;
                    Node *child = children_[k];
                    if (j != pivots[k])
                    {
                        child->data_.push_back(data_[j]);
                        child->updateRadius(dists(j, k));
                    }
                    for (unsigned int i = 0; i < degree_; ++i)
                        children_[i]->updateRange(k, dists(j, i));
                }

                for (auto &child : children_)
                {
                    // make sure degree lies between minDegree_ and maxDegree_
                    child->degree_ =
                        std::min(std::max((unsigned int)((degree_ * child->data_.size()) / data_.size()),
                                          gnat.minDegree_),
                                 gnat.maxDegree_);
                    // singleton
                    if (child->minRadius_ >= std::numeric_limits<double>::infinity())
                        child->minRadius_ = child->maxRadius_ = 0.;
#ifdef GNAT_SAMPLER
                    // set subtree size
                    child->subtreeSize_ = child->data_.size() + 1;
#endif
                }
                // this does more than clear(); it also sets capacity to 0 and frees the memory
                std::vector<_T> tmp;
                data_.swap(tmp);
                // check if new leaves need to be split
                for (auto &child : children_)
                    if (child->needToSplit(gnat))
                        child->split(gnat);
            }

            /// Insert data in nbh if it is a near neighbor. Return true iff data was added to nbh.
            bool insertNeighborK(NearQueue &nbh, std::size_t k, const _T &data, const _T &key, double dist) const
            {
                if (nbh.size() < k)
                {
                    nbh.emplace(dist, &data);
                    return true;
                }
                if (dist < nbh.top().first || (dist < std::numeric_limits<double>::epsilon() && data == key))
                {
                    nbh.pop();
                    nbh.emplace(dist, &data);
                    return true;
                }
                return false;
            }

            /// \brief Compute the k nearest neighbors of data in the tree.
            /// For k=1, isPivot is true if the nearest neighbor is a pivot
            /// (which is important during removal; removing pivots is a
            /// special case). The nodeQueue, which contains other Nodes
            /// that need to be checked for nearest neighbors, is updated.
            void nearestK(const GNAT &gnat, const _T &data, std::size_t k, NearQueue &nbh, NodeQueue &nodeQueue,
                          bool &isPivot, bool disabled = false) const
            {
                for (const auto &d : data_)
                    if (!gnat.isRemoved(d) && disabled == gnat.isDisabled(d))
                    {
                        if (insertNeighborK(nbh, k, d, data, gnat.distFun_(data, d)))
                            isPivot = false;
                    }
                if (!children_.empty())
                {
                    double dist;
                    Node *child;
                    std::size_t sz = children_.size(), offset = gnat.offset_++;
                    std::vector<double> distToPivot(sz);
                    std::vector<int> permutation(sz);
                    for (unsigned int i = 0; i < sz; ++i)
                        permutation[i] = (i + offset) % sz;

                    for (unsigned int i = 0; i < sz; ++i)
                        if (permutation[i] >= 0)
                        {
                            child = children_[permutation[i]];
                            distToPivot[permutation[i]] = gnat.distFun_(data, child->pivot_);
                            if (disabled && gnat.isDisabled(child->pivot_))
                                if (insertNeighborK(nbh, k, child->pivot_, data, distToPivot[permutation[i]]))
                                    isPivot = true;
                            if (nbh.size() == k)
                            {
                                dist = nbh.top().first;  // note difference with nearestR
                                for (unsigned int j = 0; j < sz; ++j)
                                    if (permutation[j] >= 0 && i != j &&
                                        (distToPivot[permutation[i]] - dist > child->maxRange_[permutation[j]] ||
                                         distToPivot[permutation[i]] + dist < child->minRange_[permutation[j]]))
                                        permutation[j] = -1;
                            }
                        }

                    dist = nbh.top().first;
                    for (auto p : permutation)
                        if (p >= 0)
                        {
                            child = children_[p];
                            if (nbh.size() < k || (distToPivot[p] - dist <= child->maxRadius_ &&
                                                   distToPivot[p] + dist >= child->minRadius_))
                                nodeQueue.emplace(child, distToPivot[p]);
                        }
                }
            }
            /// Insert data in nbh if it is a near neighbor.
            void insertNeighborR(NearQueue &nbh, double r, const _T &data, double dist) const
            {
                if (dist <= r)
                    nbh.emplace(dist, &data);
            }
            /// \brief Return all elements that are within distance r in nbh.
            /// The nodeQueue, which contains other Nodes that need to
            /// be checked for nearest neighbors, is updated.
            void nearestR(const GNAT &gnat, const _T &data, double r, NearQueue &nbh, NodeQueue &nodeQueue, bool disabled = false) const
            {
                double dist = r;  // note difference with nearestK

                for (const auto &d : data_)
                    if (!gnat.isRemoved(d) && disabled == gnat.isDisabled(d))
                        insertNeighborR(nbh, r, d, gnat.distFun_(data, d));
                if (!children_.empty())
                {
                    Node *child;
                    std::size_t sz = children_.size(), offset = gnat.offset_++;
                    std::vector<double> distToPivot(sz);
                    std::vector<int> permutation(sz);
                    // Not a random permutation, but processing the children in slightly different order is
                    // "good enough" to get a performance boost. A call to std::shuffle takes too long.
                    for (unsigned int i = 0; i < sz; ++i)
                        permutation[i] = (i + offset) % sz;

                    for (unsigned int i = 0; i < sz; ++i)
                        if (permutation[i] >= 0)
                        {
                            child = children_[permutation[i]];
                            distToPivot[permutation[i]] = gnat.distFun_(data, child->pivot_);
                            if (disabled && gnat.isDisabled(child->pivot_))
                                insertNeighborR(nbh, r, child->pivot_, distToPivot[permutation[i]]);
                            for (unsigned int j = 0; j < sz; ++j)
                                if (permutation[j] >= 0 && i != j &&
                                    (distToPivot[permutation[i]] - dist > child->maxRange_[permutation[j]] ||
                                     distToPivot[permutation[i]] + dist < child->minRange_[permutation[j]]))
                                    permutation[j] = -1;
                        }

                    for (auto p : permutation)
                        if (p >= 0)
                        {
                            child = children_[p];
                            if (distToPivot[p] - dist <= child->maxRadius_ &&
                                distToPivot[p] + dist >= child->minRadius_)
                                nodeQueue.emplace(child, distToPivot[p]);
                        }
                }
            }

#ifdef GNAT_SAMPLER
            double getSamplingWeight(const GNAT &gnat) const
            {
                double minR = std::numeric_limits<double>::max();
                for (auto minRange : minRange_)
                    if (minRange < minR && minRange > 0.0)
                        minR = minRange;
                minR = std::max(minR, maxRadius_);
                return std::pow(minR, gnat.estimatedDimension_) / (double)subtreeSize_;
            }
            const _T &sample(const GNAT &gnat, RNG &rng) const
            {
                if (children_.size() != 0)
                {
                    if (rng.uniform01() < 1. / (double)subtreeSize_)
                        return pivot_;
                    PDF<const Node *> distribution;
                    for (const auto &child : children_)
                        distribution.add(child, child->getSamplingWeight(gnat));
                    return distribution.sample(rng.uniform01())->sample(gnat, rng);
                }
                else
                {
                    unsigned int i = rng.uniformInt(0, data_.size());
                    return (i == data_.size()) ? pivot_ : data_[i];
                }
            }
#endif

            void list(const GNAT &gnat, std::vector<_T> &data) const
            {
                if (!gnat.isRemoved(pivot_))
                    data.push_back(pivot_);
                for (const auto &d : data_)
                    if (!gnat.isRemoved(d))
                        data.push_back(d);
                for (const auto &child : children_)
                    child->list(gnat, data);
            }

            friend std::ostream &operator<<(std::ostream &out, const Node &node)
            {
                out << "\ndegree:\t" << node.degree_;
                out << "\nminRadius:\t" << node.minRadius_;
                out << "\nmaxRadius:\t" << node.maxRadius_;
                out << "\nminRange:\t";
                for (auto minR : node.minRange_)
                    out << minR << '\t';
                out << "\nmaxRange: ";
                for (auto maxR : node.maxRange_)
                    out << maxR << '\t';
                out << "\npivot:\t" << node.pivot_;
                out << "\ndata: ";
                for (auto &data : node.data_)
                    out << data << '\t';
                out << "\nthis:\t" << &node;
#ifdef GNAT_SAMPLER
                out << "\nsubtree size:\t" << node.subtreeSize_;
                out << "\nactivity:\t" << node.activity_;
#endif
                out << "\nchildren:\n";
                for (auto &child : node.children_)
                    out << child << '\t';
                out << '\n';
                for (auto &child : node.children_)
                    out << *child << '\n';
                return out;
            }

            /// Number of child nodes
            unsigned int degree_;
            /// Data element stored in this Node
            const _T pivot_;
            /// Minimum distance between the pivot element and the elements stored in data_
            double minRadius_;
            /// Maximum distance between the pivot element and the elements stored in data_
            double maxRadius_;
            /// \brief The i-th element in minRange_ is the minimum distance between the
            /// pivot and any data_ element in the i-th child node of this node's parent.
            std::vector<double> minRange_;
            /// \brief The i-th element in maxRange_ is the maximum distance between the
            /// pivot and any data_ element in the i-th child node of this node's parent.
            std::vector<double> maxRange_;
            /// \brief The data elements stored in this node (in addition to the pivot
            /// element). An internal node has no elements stored in data_.
            std::vector<_T> data_;
            /// \brief The child nodes of this node. By definition, only internal nodes
            /// have child nodes.
            std::vector<Node *> children_;
#ifdef GNAT_SAMPLER
            /// Number of elements stored in the subtree rooted at this Node
            unsigned int subtreeSize_;
            /// \brief The extent to which a Node's maxRadius_ is increasing. A value of 0
            /// means the Node's maxRadius_ was increased the last time an element was added,
            /// while a negative value i means the Node hasn't expanded the last -i times
            /// an element was added.
            int activity_;
#endif
        };

        /// \brief The data structure containing the elements stored in this structure.
        Node *tree_{nullptr};
        /// The desired degree of each node.
        unsigned int degree_;
        /// \brief After splitting a Node, each child Node has degree equal to
        /// the default degree times the fraction of data elements from the
        /// original node that got assigned to that child Node. However, its
        /// degree can be no less than minDegree_.
        unsigned int minDegree_;
        /// \brief After splitting a Node, each child Node has degree equal to
        /// the default degree times the fraction of data elements from the
        /// original node that got assigned to that child Node. However, its
        /// degree can be no larger than maxDegree_.
        unsigned int maxDegree_;
        /// \brief Maximum number of elements allowed to be stored in a Node before
        /// it needs to be split into several nodes.
        unsigned int maxNumPtsPerLeaf_;
        /// \brief Number of elements stored in the tree.
        std::size_t size_{0};
        /// \brief If size_ exceeds rebuildSize_, the tree will be rebuilt (and
        /// automatically rebalanced), and rebuildSize_ will be doubled.
        std::size_t rebuildSize_;
        /// \brief Maximum number of removed elements that can be stored in the
        /// removed_ cache. If the cache is full, the tree will be rebuilt with
        /// the elements in removed_ actually removed from the tree.
        std::size_t removedCacheSize_;
        /// \brief The data structure used to split data into subtrees.
        GreedyKCenters<_T> pivotSelector_;
        /// \brief Cache of removed elements.
        std::unordered_set<const _T *> removed_;
        /// \brief Cache of disabled elements.
        std::unordered_set<const _T *> disabled_;
#ifdef GNAT_SAMPLER
        /// \brief Estimated dimension of the local free space.
        double estimatedDimension_;
#endif

        /// \cond IGNORE
        // used to cycle through children of a node in different orders
        mutable std::size_t offset_{0};
        /// \endcond
    };
}

#endif
