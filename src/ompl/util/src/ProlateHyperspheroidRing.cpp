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

/* Author: Shi Shenglei*/

// The class's header
#include "ompl/util/ProlateHyperspheroidRing.h"
// For OMPL exceptions
#include "ompl/util/Exception.h"
// For OMPL information
#include "ompl/util/Console.h"
// For geometric equations like prolateHyperspheroidRingMeasure
#include "ompl/util/GeometricEquations.h"

// For std::make_shared
#include <memory>

// Eigen core:
#include <Eigen/Core>

struct ompl::ProlateHyperspheroidRing::PhsrData
{
    /** \brief The measure of the PHSR. */
    double phsrMeasure_;
    /** \brief The measure of the previous ProlateHyperspheroid. */
    double pphsMeasure_;
    /** \brief Is it a ProlateHyperspheroid. */
    bool isProlateHyperspheroid;
    /** \brief The expansion factor1 */
    double expf1;
    /** \brief The expansion factor2 */
    double expf2;
    /** \brief The lalf Lebesgue measure (i.e., "volume")*/
    double expfmiddle;

    Eigen::VectorXd xCentre_;

    Eigen::VectorXd halfV;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

ompl::ProlateHyperspheroidRing::ProlateHyperspheroidRing(unsigned int n, const double focus1[], const double focus2[])
  : ProlateHyperspheroid(n, focus1, focus2), ringdataPtr_(std::make_shared<PhsrData>())
{
    // Initialize the data:
    ringdataPtr_->pphsMeasure_ = 0.0;
    ringdataPtr_->phsrMeasure_ = 0.0;
    ringdataPtr_->isProlateHyperspheroid = true;
    ringdataPtr_->expf1 = 0.0;
    ringdataPtr_->expf2 = 1.0;
    ringdataPtr_->expfmiddle = std::pow(0.5, 1.0 / static_cast<double>(n));

    Eigen::VectorXd xFocus1 = Eigen::Map<const Eigen::VectorXd>(focus1, n);
    Eigen::VectorXd xFocus2 = Eigen::Map<const Eigen::VectorXd>(focus2, n);

    ringdataPtr_->xCentre_ = 0.5 * (xFocus1 + xFocus2);
    ringdataPtr_->halfV = xFocus2 - ringdataPtr_->xCentre_;
}

void ompl::ProlateHyperspheroidRing::setTransverseDiameter(double transverseDiameter)
{
    ProlateHyperspheroid::setTransverseDiameter(transverseDiameter);

    // Store and update if changed
    ringdataPtr_->pphsMeasure_ = getPhsMeasure();
    ringdataPtr_->phsrMeasure_ = getPhsMeasure(); 
    ringdataPtr_->isProlateHyperspheroid = true;
    ringdataPtr_->expf1 = 0.0;
    ringdataPtr_->expf2 = 1.0;
    ringdataPtr_->expfmiddle = std::pow(0.5, 1.0 / static_cast<double>(getDimension()));
}

void ompl::ProlateHyperspheroidRing::setExpansionFactor(double expansionfactor)
{
    if (!isTransformUpToDate())
         throw Exception("The transverse diameter has not been set");       
    if (expansionfactor < ringdataPtr_->expf2)
        throw Exception("The expansion factor is lower than the previous one");

    unsigned int n = getDimension();
    ringdataPtr_->expf1 = ringdataPtr_->expf2;
    ringdataPtr_->expf2 = expansionfactor;
    ringdataPtr_->isProlateHyperspheroid = false;
    ringdataPtr_->expfmiddle = std::pow(0.5*( std::pow(ringdataPtr_->expf1, static_cast<double>(n)) +  
                                              std::pow(ringdataPtr_->expf2, static_cast<double>(n))),
                                        1.0 / static_cast<double>(n));

    double externalphsMeasure = std::pow(expansionfactor, static_cast<double>(n)) * getPhsMeasure();
    ringdataPtr_->phsrMeasure_ = externalphsMeasure - ringdataPtr_->pphsMeasure_;
    ringdataPtr_->pphsMeasure_ = externalphsMeasure; 
}

bool ompl::ProlateHyperspheroidRing::isInPhs(const double point[]) const
{
    if (!isTransformUpToDate())
    {
        // The transform is not up to date until the transverse diameter has been set
        throw Exception("The transverse diameter has not been set");
    }

    if (ringdataPtr_->isProlateHyperspheroid)
        return (getPathLength(point) < getTransverseDiameter());
    else
    {
         return (getPathLengthE(point, ringdataPtr_->expf1) >= ringdataPtr_->expf1 * getTransverseDiameter() && 
                 getPathLengthE(point, ringdataPtr_->expf2)  < ringdataPtr_->expf2 * getTransverseDiameter());
    }
}

bool ompl::ProlateHyperspheroidRing::isOnPhs(const double point[]) const
{
    if (!isTransformUpToDate())
    {
        // The transform is not up to date until the transverse diameter has been set
        throw Exception("The transverse diameter has not been set");
    }

     return (getPathLengthE(point, ringdataPtr_->expf2) == ringdataPtr_->expf2 * getTransverseDiameter());
}

bool ompl::ProlateHyperspheroidRing::isOutsideMiddlePhs(const double point[]) const
{
    if (!isTransformUpToDate())
    {
        // The transform is not up to date until the transverse diameter has been set
        throw Exception("The transverse diameter has not been set");
    }

     return (getPathLengthE(point, ringdataPtr_->expfmiddle) >= ringdataPtr_->expfmiddle * getTransverseDiameter());
}

double ompl::ProlateHyperspheroidRing::getPhsrMeasure() const
{
    if (!isTransformUpToDate())
    {
        // The transform is not up to date until the transverse diameter has been set, therefore we have no transverse
        // diameter and we have infinite measure
        return std::numeric_limits<double>::infinity();
    }

    return ringdataPtr_->phsrMeasure_; 
}

double ompl::ProlateHyperspheroidRing::getInternalExpansionFactor() const
{
    return ringdataPtr_->expf1;
}

double ompl::ProlateHyperspheroidRing::getExternalExpansionFactor() const
{
    return ringdataPtr_->expf2;
}

bool ompl::ProlateHyperspheroidRing::isProlateHyperspheroid() const
{
    return ringdataPtr_->isProlateHyperspheroid;
}

double ompl::ProlateHyperspheroidRing::getPathLengthE(const double point[], double expf) const
{
    if (expf == 1.0)
        return getPathLength(point);
    else 
    {
        Eigen::VectorXd x1 = ringdataPtr_->xCentre_ - expf * ringdataPtr_->halfV;
        Eigen::VectorXd x2 = ringdataPtr_->xCentre_ + expf * ringdataPtr_->halfV;

        return (x1 - Eigen::Map<const Eigen::VectorXd>(point, getDimension())).norm() +
               (Eigen::Map<const Eigen::VectorXd>(point, getDimension()) - x2).norm();
    }
}
