/*********************************************************************
 * Software License Agreement (BSD License)
 *
 *  Copyright (c) 2010, Rice University
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

/* Author: Shi Shenglei */

#ifndef OMPL_DATASTRUCTURES_GRIDNEAREST_
#define OMPL_DATASTRUCTURES_GRIDNEAREST_

#include <Eigen/Core>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <unordered_map>
#include <algorithm>
#include "ompl/util/Exception.h"
#include "ompl/datastructures/NearestNeighborsSqrtApprox.h"
#include "ompl/datastructures/NearestNeighborsGNAT.h"
#include "ompl/datastructures/NearestNeighborsGNATNoThreadSafety.h"

namespace ompl
{
    /** \brief Representation of a simple grid,
      and nearest neighbor search can be applied to _T */
    template <typename _T, typename _TT = int> // todo specialization
    class GridNearest
    {
    public:
        /// Definition of a coordinate within this grid
        using Coord = Eigen::VectorXi;

        /// Definition of a cell in this grid
        struct Cell
        {
            /// The data we store in the cell
            std::vector<_T> data;

            _TT auxData;

            /// The coordinate of the cell
            Coord coord;

            std::shared_ptr<NearestNeighbors<_T>> nn;

            Cell() = default;

            virtual ~Cell() = default;

            EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        };

        /// The datatype for arrays of cells
        using CellArray = std::vector<Cell *>;

        /** \brief The definition of a distance function */
        using DistanceFunction = std::function<double(const _T &, const _T &)>;

        /// The constructor takes the dimension of the grid as argument
        explicit GridNearest(unsigned int dimension)
        {
            setDimension(dimension);
            cnn_.reset(new NearestNeighborsGNATNoThreadSafety<Coord>()); 
            cnn_->setDistanceFunction([](const Coord &a, const Coord &b){return (a - b).template lpNorm<Eigen::Infinity>();});
        }

        /// Destructor
        virtual ~GridNearest()
        {
            freeMemory();
        }

        /// Clear all cells in the grid
        virtual void clear()
        {
            freeMemory();
            cnn_.reset();
        }

        /// Update the dimension of the grid; this should not be done
        /// unless the grid is empty
        void setDimension(unsigned int dimension)
        {
            if (!empty())
                throw;
            dimension_ = dimension;
            maxNeighbors_ = 2 * dimension_;
        }

        /// Return the dimension of the grid
        unsigned int getDimension() const
        {
            return dimension_;
        }

        void useThread(bool use)
        {
            useThread_ = use;
            if (use)
                cnn_.reset(new NearestNeighborsGNAT<Coord>()); 
            else 
                cnn_.reset(new NearestNeighborsGNATNoThreadSafety<Coord>()); 
            cnn_->setDistanceFunction([](const Coord &a, const Coord &b){return (a - b).template lpNorm<Eigen::Infinity>();});
        }

        void metricSpace(bool metric)
        {
            metric_ = metric;
        }

        /** \brief Set the distance function to use */
        virtual void setDistanceFunction(const DistanceFunction &distFun)
        {
            distFun_ = distFun;
        }

        /** \brief Get the distance function used */
        const DistanceFunction &getDistanceFunction() const
        {
            return distFun_;
        }

        /// Check if a cell exists at the specified coordinate
        bool has(const Coord &coord) const
        {
            return getCell(coord) != nullptr;
        }

        /// Get the cell at a specified coordinate
        Cell *getCell(const Coord &coord) const
        {
            auto pos = hash_.find(const_cast<Coord *>(&coord));
            Cell *c = (pos != hash_.end()) ? pos->second : nullptr;
            return c;
        }

        /** \brief Get the nearest neighbor of a point */
        virtual _T nearest(const _T &data, const Coord &coord) const
        {
            Cell *cell = getCell(coord);
            if (cell)
                return cell->nn->nearest(data);
            else if (cnn_->size() > 0) 
            {
                cell = getCell(cnn_->nearest(coord));
                return cell->nn->nearest(data);
            }
            else 
                throw Exception("No elements found in grid nearest neighbors data structure");
        }

        /** \brief Get the k-nearest neighbors of a point
         *
         * All the nearest neighbor structures currently return the neighbors in
         * sorted order, but this is not required.
         */
        virtual void nearestK(const _T &data, const Coord &coord, std::size_t k, std::vector<_T> &nbh) const
        {
            Cell *cell = getCell(coord);
            if (cell)
                return cell->nn->nearestK(data, k, nbh);
            else if (cnn_->size() > 0) 
            {
                cell = getCell(cnn_->nearest(coord));
                return cell->nn->nearestK(data, k, nbh);
            }
            else 
                throw Exception("No elements found in grid nearestK neighbors data structure");
        }

        /** \brief Get the nearest neighbors of a point, within a specified radius
         *
         * All the nearest neighbor structures currently return the neighbors in
         * sorted order, but this is not required.
         */
        virtual void nearestR(const _T &data, const Coord &coord, double radius, std::vector<_T> &nbh) const
        {
            Cell *cell = getCell(coord);
            if (cell)
                return cell->nn->nearestR(data, radius, nbh);
            else if (cnn_->size() > 0) 
            {
                cell = getCell(cnn_->nearest(coord));
                return cell->nn->nearestR(data, radius, nbh);
            }
            else 
                throw Exception("No elements found in grid nearestR neighbors data structure");
        }

        /// Get the list of neighbors for a given cell
        void neighbors(const Cell *cell, CellArray &list) const
        {
            Coord test = cell->coord;
            neighbors(test, list);
        }

        /// Get the list of neighbors for a given coordinate
        void neighbors(const Coord &coord, CellArray &list) const
        {
            Coord test = coord;
            neighbors(test, list);
        }

        /// Get the list of neighbors for a given coordinate
        void neighbors(Coord &coord, CellArray &list) const
        {
            list.reserve(list.size() + maxNeighbors_);

            for (int i = dimension_ - 1; i >= 0; --i)
            {
                coord[i]--;

                auto pos = hash_.find(&coord);
                Cell *cell = (pos != hash_.end()) ? pos->second : nullptr;

                if (cell)
                    list.push_back(cell);
                coord[i] += 2;

                pos = hash_.find(&coord);
                cell = (pos != hash_.end()) ? pos->second : nullptr;

                if (cell)
                    list.push_back(cell);
                coord[i]--;
            }
        }

        /// Get the connected components formed by the cells in this grid (based on neighboring relation)
        std::vector<std::vector<Cell *>> components() const
        {
            using ComponentHash = std::unordered_map<Coord *, int, HashFunCoordPtr, EqualCoordPtr>;

            int components = 0;
            ComponentHash ch;
            std::vector<std::vector<Cell *>> res;

            for (auto & i: hash_)
            {
                Cell *c0 = i.second;
                auto pos = ch.find(&c0->coord);
                int comp = (pos != ch.end()) ? pos->second : -1;

                if (comp < 0)
                {
                    res.resize(res.size() + 1);
                    std::vector<Cell *> &q = res.back();
                    q.push_back(c0);
                    std::size_t index = 0;
                    while (index < q.size())
                    {
                        Cell *c = q[index++];
                        pos = ch.find(&c->coord);
                        comp = (pos != ch.end()) ? pos->second : -1;

                        if (comp < 0)
                        {
                            ch.insert(std::make_pair(&c->coord, components));
                            std::vector<Cell *> nbh;
                            neighbors(c, nbh);
                            for (const auto &n : nbh)
                            {
                                pos = ch.find(&n->coord);
                                comp = (pos != ch.end()) ? pos->second : -1;
                                if (comp < 0)
                                    q.push_back(n);
                            }
                        }
                        else
                        {
                            --index;
                            q.erase(q.begin() + index);
                        }
                    }
                    ++components;
                }
            }
            std::sort(res.begin(), res.end(), SortComponents());
            return res;
        }

        /// Instantiate a new cell at given coordinates; optionally
        /// Return the list of future neighbors.
        /// Note: this call only creates the cell, but does not add it to the grid.
        /// It however updates the neighbor count for neighboring cells
        virtual Cell *createCell(const Coord &coord, CellArray *nbh = nullptr)
        {
            Cell *cell = new Cell();
            cell->coord = coord;
            if (metric_)
            {
                if (useThread_)
                    cell->nn.reset(new NearestNeighborsGNAT<_T>());
                else
                    cell->nn.reset(new NearestNeighborsGNATNoThreadSafety<_T>());
            }
            else 
                cell->nn.reset(new NearestNeighborsSqrtApprox<_T>());
            cell->nn->setDistanceFunction(distFun_);
            if (nbh)
                neighbors(cell->coord, *nbh);
            return cell;
        }

        /// Remove a cell from the grid. If the cell has not been
        /// Added to the grid, only update the neighbor list
        virtual bool remove(Cell *cell)
        {
            if (cell)
            {
                auto pos = hash_.find(&cell->coord);
                if (pos != hash_.end())
                {
                    hash_.erase(pos);
                    cnn_->remove(cell->coord);
                    return true;
                }
            }
            return false;
        }

        /// Add an instantiated cell to the grid
        virtual void add(Cell *cell)
        {
            hash_.insert(std::make_pair(&cell->coord, cell));
            cnn_->add(cell->coord);
        }

        virtual void add(const _T &data, const Coord &coord)
        {
            Cell *cell = getCell(coord);
            if (!cell)
            {
                cell = createCell(coord);
                add(cell);
            }
            add(data, cell);
        }

        virtual void add(const _T &data, Cell *cell)
        {
            cell->data.push_back(data);
            cell->nn->add(data);
        }

        virtual bool remove(const _T &data, const Coord &coord)
        {
            Cell *cell = getCell(coord);
            if (cell)
                return remove(data, cell);
            else 
                throw Exception("The to be removed element is not in the grid nearest neighbors data structure");
            return false;
        }

        virtual bool remove(const _T &data, Cell *cell)
        {
            bool found = false;
            for (std::size_t i = cell->data.size() - 1; i < cell->data.size(); i--)
            {
                if (cell->data[i] == data)
                {
                    std::iter_swap(cell->data.begin() + i, cell->data.end() - 1);
                    cell->data.pop_back();
                    cell->nn->remove(data);
                    found = true;
                    break;
                }
            }
            return found;
        }

        /// Clear the memory occupied by a cell; do not call this function unless remove() was called first
        virtual void destroyCell(Cell *cell) const
        {
            delete cell;
        }

        /// Get the data stored in the cells we are aware of
        void getContent(std::vector<std::vector<_T>> &content) const
        {
            for (const auto &h : hash_)
                content.push_back(h.second->data);
        }

        /// Get the set of coordinates where there are cells
        void getCoordinates(std::vector<Coord *> &coords) const
        {
            for (const auto &h : hash_)
                coords.push_back(h.first);
        }

        /// Get the set of instantiated cells in the grid
        void getCells(CellArray &cells) const
        {
            for (const auto &h : hash_)
                cells.push_back(h.second);
        }

        /// Print the value of a coordinate to a stream
        void printCoord(Coord &coord, std::ostream &out = std::cout) const
        {
            out << "[ ";
            for (unsigned int i = 0; i < dimension_; ++i)
                out << coord[i] << " ";
            out << "]" << std::endl;
        }

        /// Check if the grid is empty
        bool empty() const
        {
            return hash_.empty();
        }

        /// Check the size of the grid
        unsigned int size() const
        {
            return hash_.size();
        }

        /// Print information about the data in this grid structure
        virtual void status(std::ostream &out = std::cout) const
        {
            out << size() << " total cells " << std::endl;
            const std::vector<std::vector<Cell *>> &comp = components();
            out << comp.size() << " connected components: ";
            for (const auto &c : comp)
                out << c.size() << " ";
            out << std::endl;
        }

    protected:
        /// Free the allocated memory
        void freeMemory()
        {
            CellArray content;
            getCells(content);
            hash_.clear();
            for (auto &c : content)
                delete c;
        }

        /// Hash function for coordinates; see
        /// http://www.cs.hmc.edu/~geoff/classes/hmc.cs070.200101/homework10/hashfuncs.html
        struct HashFunCoordPtr
        {
            /// Hash function for coordinates
            std::size_t operator()(const Coord *const s) const
            {
                unsigned long h = 0;
                for (int i = s->size() - 1; i >= 0; --i)
                {
                    int high = h & 0xf8000000;
                    h = h << 5;
                    h = h ^ (high >> 27);
                    h = h ^ std::abs((*s)[i]);
                }
                return (std::size_t)h;
            }
        };

        /// Equality operator for coordinate pointers
        struct EqualCoordPtr
        {
            /// Equality operator for coordinate pointers
            bool operator()(const Coord *const c1, const Coord *const c2) const
            {
                return *c1 == *c2;
            }
        };

        /// Define the datatype for the used hash structure
        using CoordHash = std::unordered_map<Coord *, Cell *, HashFunCoordPtr, EqualCoordPtr>;

        /// Helper to sort components by size
        struct SortComponents
        {
            /// Helper to sort components by size
            bool operator()(const std::vector<Cell *> &a, const std::vector<Cell *> &b) const
            {
                return a.size() > b.size();
            }
        };

    public:
        /// We only allow const iterators
        using iterator = typename CoordHash::const_iterator;

        /// Return the begin() iterator for the grid
        iterator begin() const
        {
            return hash_.begin();
        }

        /// Return the end() iterator for the grid
        iterator end() const
        {
            return hash_.end();
        }

        iterator erase(iterator it)
        {
            it = hash_.erase(it);
            return it;
        }

    protected:
        /// The dimension of the grid
        unsigned int dimension_;

        /// The maximum number of neighbors a cell can have (2 * dimension)
        unsigned int maxNeighbors_;

        /// The hash holding the cells
        CoordHash hash_;

        std::shared_ptr<NearestNeighbors<Coord>> cnn_;

    protected:
        /** \brief The used distance function */
        DistanceFunction distFun_;

        bool useThread_{false};

        bool metric_{true};
    };
}  // namespace ompl

#endif
