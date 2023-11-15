/*********************************************************************
* Software License Agreement (BSD License)
*
*  Copyright (c) 2008, Willow Garage, Inc.
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
*   * Neither the name of the Willow Garage nor the names of its
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

/* Author: Ioan Sucan, Shi Shenglei */

#ifndef OMPL_DATASTRUCTURES_NEAREST_NEIGHBORS_LINEAR_
#define OMPL_DATASTRUCTURES_NEAREST_NEIGHBORS_LINEAR_

#include "ompl/datastructures/NearestNeighbors.h"
#include "ompl/util/Exception.h"
#include <algorithm>

namespace ompl
{
    /** \brief A nearest neighbors datastructure that uses linear
        search.

        \li Search for nearest neighbor is O(n).
        \li Search for k-nearest neighbors is  O(n log(k)).
        \li Search for neighbors within a range is O(n log(n)).
        \li Adding an element to the datastructure is O(1).
        \li Removing an element from the datastructure O(n).
    */
    template <typename _T>
    class NearestNeighborsLinear : public NearestNeighbors<_T>
    {
    public:
        NearestNeighborsLinear() : NearestNeighbors<_T>()
        {
        }

        ~NearestNeighborsLinear() override = default;

        void clear() override
        {
            data_.clear();
            disabled_ = 0;
        }

        bool reportsSortedResults() const override
        {
            return true;
        }

        void add(const _T &data) override
        {
            data_.push_back(data);
            if (disabled_)
                std::iter_swap(data_.end() - 1, data_.end() - 1 - disabled_);
        }

        void add(const std::vector<_T> &data) override
        {
            data_.reserve(data_.size() + data.size());
            if (disabled_)
            {
                for (auto & d : data)
                {
                    data_.push_back(d);
                    std::iter_swap(data_.end() - 1, data_.end() - 1 - disabled_);
                }
            }
            else 
                data_.insert(data_.end(), data.begin(), data.end());
        }

        void enable(const _T &data) override
        {
            for (auto it = data_.end() - disabled_; it != data_.end(); it++)
            {
                if ((*it) == data)
                {
                    std::iter_swap(it, data_.end() - disabled_);
                    disabled_--;
                    break;
                }
            }
        }

        void disable(const _T &data) override
        {
            for (std::size_t i = data_.size() - disabled_ - 1; i < data_.size(); i--)
            {
                if (data_[i] == data)
                {
                    disabled_++;
                    std::iter_swap(data_.begin() + i, data_.end() - disabled_);
                    break;
                }
            }
        }

        bool remove(const _T &data) override
        {
            for (std::size_t i = data_.size() - 1; i < data_.size(); i--)
                if (data_[i] == data)
                {
                    if (!disabled_ || i >= data_.size() - disabled_)
                    {
                        std::iter_swap(data_.begin() + i, data_.end() - 1);
                        data_.pop_back();
                        if (disabled_)
                            disabled_--;
                    }
                    else 
                        data_.erase(data_.begin() + i);
                    return true;
                }
            return false;
        }

        _T nearest(const _T &data) const override
        {
            const std::size_t sz = data_.size() - disabled_;
            std::size_t pos = sz;
            double dmin = 0.0;
            for (std::size_t i = 0; i < sz; ++i)
            {
                double distance = NearestNeighbors<_T>::distFun_(data_[i], data);
                if (pos == sz || dmin > distance)
                {
                    pos = i;
                    dmin = distance;
                }
            }
            if (pos != sz)
                return data_[pos];

            throw Exception("No elements found in nearest neighbors data structure");
        }

        /// Return the k nearest neighbors in sorted order
        void nearestK(const _T &data, std::size_t k, std::vector<_T> &nbh) const override
        {
            nbh = data_;
            nbh.resize(data_.size() - disabled_);
            if (nbh.size() > k)
            {
                std::partial_sort(nbh.begin(), nbh.begin() + k, nbh.end(),
                                  ElemSort(data, NearestNeighbors<_T>::distFun_));
                nbh.resize(k);
            }
            else
            {
                std::sort(nbh.begin(), nbh.end(), ElemSort(data, NearestNeighbors<_T>::distFun_));
            }
        }

        /// Return the nearest neighbors within distance \c radius in sorted order
        void nearestR(const _T &data, double radius, std::vector<_T> &nbh) const override
        {
            nbh.clear();
            for (std::size_t i = 0; i < data_.size() - disabled_; i++)
            {
                _T d = data_[i];
                if (NearestNeighbors<_T>::distFun_(d, data) <= radius)
                    nbh.push_back(d);
            }
            std::sort(nbh.begin(), nbh.end(), ElemSort(data, NearestNeighbors<_T>::distFun_));
        }

        // Get the nearest neighbor of a point in the disabled set */
        _T nearestDisabled(const _T &data) const override
        {
            const std::size_t sz = data_.size();
            std::size_t pos = sz;
            double dmin = 0.0;
            for (std::size_t i = sz - disabled_; i < sz; ++i)
            {
                double distance = NearestNeighbors<_T>::distFun_(data_[i], data);
                if (pos == sz || dmin > distance)
                {
                    pos = i;
                    dmin = distance;
                }
            }
            if (pos != sz)
                return data_[pos];

            throw Exception("No elements found in disabled nearest neighbors data structure");
        }

        // Get the k-nearest neighbors of a point in the disabled set
        void nearestKDisabled(const _T &data, std::size_t k, std::vector<_T> &nbh) const override 
        {
            nbh.clear();
            nbh.insert(nbh.end(), data_.end() - disabled_, data_.end());
            if (nbh.size() > k)
            {
                std::partial_sort(nbh.begin(), nbh.begin() + k, nbh.end(),
                                  ElemSort(data, NearestNeighbors<_T>::distFun_));
                nbh.resize(k);
            }
            else
            {
                std::sort(nbh.begin(), nbh.end(), ElemSort(data, NearestNeighbors<_T>::distFun_));
            }
        }

        // Get the nearest neighbors of a point in the disabled set, within a specified radius
        void nearestRDisabled(const _T &data, double radius, std::vector<_T> &nbh) const override
        {
            nbh.clear();
            for (std::size_t i = data_.size() - disabled_; i < data_.size(); i++)
            {
                _T d = data_[i];
                if (NearestNeighbors<_T>::distFun_(d, data) <= radius)
                    nbh.push_back(d);
            }
            std::sort(nbh.begin(), nbh.end(), ElemSort(data, NearestNeighbors<_T>::distFun_));
        }

        std::size_t size() const override
        {
            return data_.size();
        }

        std::size_t disabledSize() const override
        {
            return disabled_;
        }

        void list(std::vector<_T> &data) const override
        {
            data = data_;
        }

    protected:
        /** \brief The data elements stored in this structure */
        std::vector<_T> data_;
        /** \brief The disabled data elements in this structure, they are always at the end of data_ */
        std::size_t disabled_{0};

    private:
        struct ElemSort
        {
            ElemSort(const _T &e, const typename NearestNeighbors<_T>::DistanceFunction &df) : e_(e), df_(df)
            {
            }

            bool operator()(const _T &a, const _T &b) const
            {
                return df_(a, e_) < df_(b, e_);
            }

            const _T &e_;
            const typename NearestNeighbors<_T>::DistanceFunction &df_;
        };
    };
}

#endif
