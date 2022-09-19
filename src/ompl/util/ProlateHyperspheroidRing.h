/*********************************************************************
* Software License Agreement (BSD License)
*
*  Copyright (c) 2014, University of Toronto
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
*   * Neither the name of the University of Toronto nor the names of its
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

#ifndef OMPL_UTIL_PROLATE_HYPERSPHEROIDRING_
#define OMPL_UTIL_PROLATE_HYPERSPHEROIDRING_

#include <ompl/util/ProlateHyperspheroid.h>

namespace ompl
{
    // A forward declaration of the prolate hyperspheroid ring class
    OMPL_CLASS_FORWARD(ProlateHyperspheroidRing);

    class ProlateHyperspheroidRing : public ProlateHyperspheroid
    {
    public:
        /** \brief The description of an n-dimensional prolate hyperspheroid ring */
        ProlateHyperspheroidRing(unsigned int n, const double focus1[], const double focus2[]);

        void setTransverseDiameter(double transverseDiameter);

        /** \brief Set the expansion factor of the PHSr */
        void setExpansionFactor(double expansionfactor);

        /** \brief Check if the given point lies \e in the PHSR. */
        bool isInPhs(const double point[]) const;

        /** \brief Check if the given point lies \e on the PHSR. */
        bool isOnPhs(const double point[]) const;

        /** \brief Check if the given point lies outside the half volume of the PHSR. */
        bool isOutsideMiddlePhs(const double point[]) const;

        /** \brief The measure of the PHSR */
        double getPhsrMeasure() const;

        double getInternalExpansionFactor() const;
        double getExternalExpansionFactor() const;
        bool   isProlateHyperspheroid() const;

        double getPathLengthE(const double point[], double expf) const;

    protected:
    private:
        /** \brief A forward declaration to the data structure class for the phsr. */
        struct PhsrData;

        /** \brief A shared pointer to the actual data of a ProlateHyperspheroidRing. Used to hide Eigen from the header. */
        std::shared_ptr<PhsrData> ringdataPtr_;
    };
}

#endif
