/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, development version     *
*                (c) 2006-2019 INRIA, USTL, UJF, CNRS, MGH                    *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This program is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this program. If not, see <http://www.gnu.org/licenses/>.        *
*******************************************************************************
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#pragma once

#include <sofa/core/collision/Intersection.h>
#include <SofaImplicitField/components/collision/ScalarFieldCollisionModel.h>
#include <SofaBaseCollision/DiscreteIntersection.h>
#include <SofaMeshCollision/PointModel.h>
#include <SofaMeshCollision/LineModel.h> 

namespace sofa::component::collision
{
namespace _scalarfieldintersectionline_
{
using sofa::core::collision::BaseIntersector;
using sofa::component::collision::ScalarFieldCollisionElement;
using sofa::component::collision::DiscreteIntersection;
using sofa::component::collision::PointModel;
using sofa::component::collision::LineModel;
using sofa::defaulttype::Vec3d ;


class SOFA_SOFAIMPLICITFIELD_API ScalarFieldIntersectionLine : public BaseIntersector
{
    typedef DiscreteIntersection::OutputVector OutputVector;

public:
    ScalarFieldIntersectionLine(DiscreteIntersection* object);

    bool testIntersection(ScalarFieldCollisionElement&, Line&);
    int computeIntersection(ScalarFieldCollisionElement&, Line&, OutputVector*);
    ///double Interpolation(double a, double b , double lambda);
    Vec3d Interpolation_point(Vec3d& P1 , Vec3d& P2 , double lambda);
    int n;


protected:
    DiscreteIntersection* intersection;
};
} /// namespace _scalarfieldintersectionline_
} /// sofa::component::collision
