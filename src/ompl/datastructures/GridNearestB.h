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

/* Author: Shi Shenglei */

#ifndef OMPL_DATASTRUCTURES_GRIDNEAREST_B_
#define OMPL_DATASTRUCTURES_GRIDNEAREST_B_

#include "ompl/datastructures/GridNearestN.h"
#include "ompl/datastructures/BinaryHeap.h"
#include "ompl/util/DisableCompilerWarning.h"

OMPL_PUSH_DISABLE_CLANG_WARNING(-Woverloaded-virtual)

namespace ompl
{
    /** \brief This class defines a grid that keeps track of its boundary:
     * it distinguishes between interior and exterior cells.  */
    template <typename _T, typename _TT, class LessThanExternal = std::less<_TT>, class LessThanInternal = LessThanExternal>
    class GridNearestB : public GridNearestN<_T, _TT>
    {
    public:
        /// Definition of a cell in this grid
        using Cell = typename GridNearestN<_T, _TT>::Cell;

        /// The datatype for arrays of cells
        using CellArray = typename GridNearestN<_T, _TT>::CellArray;

        /// Datatype for cell coordinates
        using Coord = typename GridNearestN<_T, _TT>::Coord;

    protected:
        /// \cond IGNORE
        // the type of cell here needs an extra pointer to allow the updatable heap to work fast
        // however, this stays hidden from the user
        struct CellX : public Cell
        {
            CellX() : Cell()
            {
            }

            ~CellX() override = default;

            void *heapElement;

            EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        };

        /// \endcond

    public:
        /// Event to be called when a cell's priority is to be updated
        using EventCellUpdate = void (*)(Cell *, void *);

        /// Constructor
        explicit GridNearestB(unsigned int dimension) : GridNearestN<_T, _TT>(dimension)
        {
            setupHeaps();
        }

        ~GridNearestB() override
        {
            clearHeaps();
        }

        /// Set the function callback and to be called when a cell's
        /// priority is updated
        void onCellUpdate(EventCellUpdate event, void *arg)
        {
            eventCellUpdate_ = event;
            eventCellUpdateData_ = arg;
        }

        /// Return the cell that is at the top of the heap maintaining internal cells
        Cell *topInternal() const
        {
            auto *top = static_cast<Cell *>(internal_.top()->data);
            return top ? top : topExternal();
        }

        /// Return the cell that is at the top of the heap maintaining external cells
        Cell *topExternal() const
        {
            auto *top = static_cast<Cell *>(external_.top()->data);
            return top ? top : topInternal();
        }

        void getHeapContent(CellArray &cells) const 
        {
            std::vector<CellX *> contentsI, contentsE;
            internal_.getContent(contentsI);
            external_.getContent(contentsE);

            cells.clear();
            cells.reserve(contentsI.size() + contentsE.size());

            for (auto & c : contentsI)
                cells.push_back(static_cast<Cell *>(c));
            for (auto & c : contentsE)
                cells.push_back(static_cast<Cell *>(c));
        }

        /// Return the number of internal cells
        unsigned int countInternal() const
        {
            return internal_.size();
        }

        /// Return the number of external cells
        unsigned int countExternal() const
        {
            return external_.size();
        }

        /// Return the fraction of external cells
        double fracExternal() const
        {
            return external_.empty() ? 0.0 : (double)(external_.size()) / (double)(external_.size() + internal_.size());
        }

        /// Return the fraction of internal cells
        double fracInternal() const
        {
            return 1.0 - fracExternal();
        }

        /// Update the position in the heaps for a particular cell.
        void update(Cell *cell)
        {
            eventCellUpdate_(cell, eventCellUpdateData_);
            if (cell->border)
                external_.update(
                    reinterpret_cast<typename externalBHeap::Element *>(static_cast<CellX *>(cell)->heapElement));
            else
                internal_.update(
                    reinterpret_cast<typename internalBHeap::Element *>(static_cast<CellX *>(cell)->heapElement));
        }

        /// Update all cells and reconstruct the heaps
        void updateAll()
        {
            std::vector<Cell *> cells;
            this->getCells(cells);
            for (int i = cells.size() - 1; i >= 0; --i)
                eventCellUpdate_(cells[i], eventCellUpdateData_);
            external_.rebuild();
            internal_.rebuild();
        }

        /// Create a cell but do not add it to the grid; update neighboring cells however
        virtual Cell *createCell(const Coord &coord, CellArray *nbh = nullptr)
        {
            auto *cell = new CellX();
            cell->coord = coord;
            if (GridNearestN<_T, _TT>::metric_)
            {
                if (GridNearestN<_T, _TT>::useThread_)
                    cell->nn.reset(new NearestNeighborsGNAT<_T>());
                else
                    cell->nn.reset(new NearestNeighborsGNATNoThreadSafety<_T>());
            }
            else 
                cell->nn.reset(new NearestNeighborsSqrtApprox<_T>());
            cell->nn->setDistanceFunction(GridNearestN<_T, _TT>::distFun_);

            CellArray *list = nbh ? nbh : new CellArray();
            this->neighbors(cell->coord, *list);

            std::vector<bool> wasBorder;
            cell->nnbh.reserve(list->size());
            wasBorder.reserve(list->size());
            for (auto & cl : *list)
            {
                cell->nnbh.push_back(cl);
                cl->nnbh.push_back(cell);
                wasBorder.push_back(cl->border);
                cl->neighbors++;
                if (cl->border && cl->neighbors >= GridNearestN<_T, _TT>::interiorCellNeighborsLimit_)
                    cl->border = false;
            }

            if (GridNearestN<_T, _TT>::removed_)
            {
                for (auto & cl : *list)
                {
                    std::size_t i =0;
                    if (!cl->removed)
                    {
                        auto *c = static_cast<CellX *>(cl);
                        eventCellUpdate_(c, eventCellUpdateData_);
                        if (c->border)
                            external_.update(reinterpret_cast<typename externalBHeap::Element *>(c->heapElement));
                        else
                        {
                            if (wasBorder[i])
                            {
                                external_.remove(reinterpret_cast<typename externalBHeap::Element *>(c->heapElement));
                                internal_.insert(c);
                            }
                            else
                                internal_.update(reinterpret_cast<typename internalBHeap::Element *>(c->heapElement));
                        }
                    }
                    i++;
                }
            }
            else 
            {
                for (auto & cl : *list)
                {
                    std::size_t i =0;
                    auto *c = static_cast<CellX *>(cl);
                    eventCellUpdate_(c, eventCellUpdateData_);
                    if (c->border)
                        external_.update(reinterpret_cast<typename externalBHeap::Element *>(c->heapElement));
                    else
                    {
                        if (wasBorder[i])
                        {
                            external_.remove(reinterpret_cast<typename externalBHeap::Element *>(c->heapElement));
                            internal_.insert(c);
                        }
                        else
                            internal_.update(reinterpret_cast<typename internalBHeap::Element *>(c->heapElement));
                    }
                    i++;
                }
            }

            cell->neighbors = GridNearestN<_T, _TT>::numberOfBoundaryDimensions(cell->coord);
            if (GridNearestN<_T, _TT>::removed_)
                cell->neighbors += std::accumulate(list->begin(), list->end(), 0, [](int a, Cell *b){ int bb = b->removed ? 0 : 1; return a + bb;});
            else
                cell->neighbors += list->size();
            if (cell->border && cell->neighbors >= GridNearestN<_T, _TT>::interiorCellNeighborsLimit_)
                cell->border = false;

            if (!nbh)
                delete list;

            return static_cast<Cell *>(cell);
        }

        /// Add the cell to the grid
        virtual void add(Cell *cell)
        {
            auto *ccell = static_cast<CellX *>(cell);
            eventCellUpdate_(ccell, eventCellUpdateData_);

            GridNearestN<_T, _TT>::add(cell);

            if (cell->border)
                external_.insert(ccell);
            else
                internal_.insert(ccell);
        }

        /// Remove a cell from the grid
        virtual bool remove(Cell *cell)
        {
            if (!cell)
                return false;
            std::vector<bool> wasBorder;
            wasBorder.reserve(cell->nnbh.size());
            for (auto & cl : cell->nnbh)
            {
                wasBorder.push_back(cl->border);
                cl->neighbors--;
                if (!cl->border && cl->neighbors < GridNearestN<_T, _TT>::interiorCellNeighborsLimit_)
                    cl->border = true;
                for (std::size_t i = cl->nnbh.size() - 1; i < cl->nnbh.size(); i--)
                {
                    if (cl->nnbh[i] == cell)
                    {
                        std::iter_swap(cl->nnbh.begin() + i, cl->nnbh.end() - 1);
                        cl->nnbh.pop_back();
                        break;
                    }
                }
            }

            if (GridNearestN<_T, _TT>::removed_)
            {
                for (auto & cl : cell->nnbh)
                {
                    std::size_t i = 0;
                    if (!cl->removed)
                    {
                        CellX *c = static_cast<CellX *>(cl);
                        eventCellUpdate_(c, eventCellUpdateData_);
                        if (c->border)
                        {
                            if (wasBorder[i])
                                external_.update(reinterpret_cast<typename externalBHeap::Element *>(c->heapElement));
                            else
                            {
                                internal_.remove(reinterpret_cast<typename internalBHeap::Element *>(c->heapElement));
                                external_.insert(c);
                            }
                        }
                        else
                            internal_.update(reinterpret_cast<typename internalBHeap::Element *>(c->heapElement));
                    }
                    i++;
                }
            }
            else 
            {
                for (auto & cl : cell->nnbh)
                {
                    std::size_t i =0;
                    CellX *c = static_cast<CellX *>(cl);
                    eventCellUpdate_(c, eventCellUpdateData_);
                    if (c->border)
                    {
                        if (wasBorder[i])
                            external_.update(reinterpret_cast<typename externalBHeap::Element *>(c->heapElement));
                        else
                        {
                            internal_.remove(reinterpret_cast<typename internalBHeap::Element *>(c->heapElement));
                            external_.insert(c);
                        }
                    }
                    else
                        internal_.update(reinterpret_cast<typename internalBHeap::Element *>(c->heapElement));
                    i++;
                }
            }

            for (std::size_t i = 1; i < cell->nbh.size(); i++)
            {
                for (std::size_t j = 0; j < cell->nbh[i].size(); j++)
                {
                    Cell* cl = cell->nbh[i][j];
                    for (std::size_t k = 0; k < cl->nbh[i].size(); k++)
                    {
                        if (cl->nbh[i][k] == cell)
                        {
                            cl->nbh[i].erase(cl->nbh[i].begin() + k);
                            break;
                        }
                    }
                }
            }

            auto pos = GridNearestN<_T, _TT>::hash_.find(&cell->coord);
            if (pos != GridNearestN<_T, _TT>::hash_.end())
            {
                GridNearestN<_T, _TT>::hash_.erase(pos);
                GridNearestN<_T, _TT>::cnn_->remove(cell->coord);
                auto *cx = static_cast<CellX *>(cell);
                if (cx->border)
                    external_.remove(reinterpret_cast<typename externalBHeap::Element *>(cx->heapElement));
                else
                    internal_.remove(reinterpret_cast<typename internalBHeap::Element *>(cx->heapElement));
                return true;
            }
            return false;
        }

        virtual void add(const _T &data, const Coord &coord)
        {
            Cell *cell = GridNearestN<_T, _TT>::getCell(coord);
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
            Cell *cell = GridNearestN<_T, _TT>::getCell(coord);
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

        /// Remove a cell from the grid
        virtual bool lazyRemove(Cell *cell)
        {
            if (!cell)
                return false;
            if (cell->removed)
                return false;
            for (auto & cl : cell->nnbh)
            {
                bool wasBorder = cl->border;
                cl->neighbors--;
                if (!cl->border && cl->neighbors < GridNearestN<_T, _TT>::interiorCellNeighborsLimit_)
                    cl->border = true;

                if (!cl->removed)
                {
                    CellX *c = static_cast<CellX *>(cl);
                    eventCellUpdate_(c, eventCellUpdateData_);
                    if (c->border)
                    {
                        if (wasBorder)
                            external_.update(reinterpret_cast<typename externalBHeap::Element *>(c->heapElement));
                        else
                        {
                            internal_.remove(reinterpret_cast<typename internalBHeap::Element *>(c->heapElement));
                            external_.insert(c);
                        }
                    }
                    else
                        internal_.update(reinterpret_cast<typename internalBHeap::Element *>(c->heapElement));
                }
            }

            auto pos = GridNearestN<_T, _TT>::hash_.find(&cell->coord);
            if (pos != GridNearestN<_T, _TT>::hash_.end())
            {
                cell->removed = true;
                GridNearestN<_T, _TT>::removed_++;
                GridNearestN<_T, _TT>::cnn_->remove(cell->coord);
                auto *cx = static_cast<CellX *>(cell);
                if (cx->border)
                    external_.remove(reinterpret_cast<typename externalBHeap::Element *>(cx->heapElement));
                else
                    internal_.remove(reinterpret_cast<typename internalBHeap::Element *>(cx->heapElement));
                return true;
            }
            return false;
        }

        /// Lazily add the cell to the grid
        virtual void lazyAdd(Cell *cell)
        {
            if (!cell)
                return;
            if (!cell->removed)
                return;
            for (auto & cl : cell->nnbh)
            {
                bool wasBorder = cl->border;
                cl->neighbors++;
                if (cl->border && cl->neighbors >= GridNearestN<_T, _TT>::interiorCellNeighborsLimit_)
                    cl->border = false;

                if (!cl->removed)
                {
                    CellX *c = static_cast<CellX *>(cl);
                    eventCellUpdate_(c, eventCellUpdateData_);
                    if (c->border)
                        external_.update(reinterpret_cast<typename externalBHeap::Element *>(c->heapElement));
                    else
                    {
                        if (wasBorder)
                        {
                            external_.remove(reinterpret_cast<typename externalBHeap::Element *>(c->heapElement));
                            internal_.insert(c);
                        }
                        else
                            internal_.update(reinterpret_cast<typename internalBHeap::Element *>(c->heapElement));
                    }
                }
            }

            auto *ccell = static_cast<CellX *>(cell);
            eventCellUpdate_(ccell, eventCellUpdateData_);
            if (cell->border)
                external_.insert(ccell);
            else
                internal_.insert(ccell);

            auto pos = GridNearestN<_T, _TT>::hash_.find(&cell->coord);
            if (pos != GridNearestN<_T, _TT>::hash_.end())
            {
                cell->removed = false;
                GridNearestN<_T, _TT>::removed_--;
                GridNearestN<_T, _TT>::cnn_->add(cell->coord);
            }
        }

        void clear() override
        {
            GridNearestN<_T, _TT>::clear();
            clearHeaps();
        }

        void status(std::ostream &out = std::cout) const override
        {
            GridNearestN<_T, _TT>::status(out);
            out << countInternal() << " internal cells" << std::endl;
            out << countExternal() << " external cells" << std::endl;
        }

    protected:
        /// Pointer to function to be called when a cell needs to be updated
        EventCellUpdate eventCellUpdate_;

        /// Data to be passed to function pointer above
        void *eventCellUpdateData_;

        /// Default no-op update routine for a cell
        static void noCellUpdate(Cell * /*unused*/, void * /*unused*/)
        {
        }

        /// Set the update procedure for the heaps of internal and external cells
        void setupHeaps()
        {
            eventCellUpdate_ = &noCellUpdate;
            eventCellUpdateData_ = nullptr;
            internal_.onAfterInsert(&setHeapElementI, nullptr);
            external_.onAfterInsert(&setHeapElementE, nullptr);
        }

        /// Clear the data from both heaps
        void clearHeaps()
        {
            internal_.clear();
            external_.clear();
        }

        /// Define order for internal cells
        struct LessThanInternalCell
        {
            bool operator()(const CellX *const a, const CellX *const b) const
            {
                return lt_(a->auxData, b->auxData);
            }

        private:
            LessThanInternal lt_;
        };

        /// Define order for external cells
        struct LessThanExternalCell
        {
            bool operator()(const CellX *const a, const CellX *const b) const
            {
                return lt_(a->auxData, b->auxData);
            }

        private:
            LessThanExternal lt_;
        };

        /// Datatype for a heap of cells containing interior cells
        using internalBHeap = BinaryHeap<CellX *, LessThanInternalCell>;

        /// Datatype for a heap of cells containing exterior cells
        using externalBHeap = BinaryHeap<CellX *, LessThanExternalCell>;

        /// Routine used internally for keeping track of binary heap elements for internal cells
        static void setHeapElementI(typename internalBHeap::Element *element, void * /*unused*/)
        {
            element->data->heapElement = reinterpret_cast<void *>(element);
        }

        /// Routine used internally for keeping track of binary heap elements for external cells
        static void setHeapElementE(typename externalBHeap::Element *element, void * /*unused*/)
        {
            element->data->heapElement = reinterpret_cast<void *>(element);
        }

        /// The heap of interior cells
        internalBHeap internal_;

        /// The heap of external cells
        externalBHeap external_;
    };
}

OMPL_POP_CLANG

#endif
