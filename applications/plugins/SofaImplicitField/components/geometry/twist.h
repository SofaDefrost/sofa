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
#ifndef SOFA_IMPLICIT_CUBEFIELD_H
#define SOFA_IMPLICIT_CUBEFIELD_H

#include "ScalarField.h"

namespace sofa
{

namespace component
{

namespace geometry
{

namespace _twist_
{

using sofa::defaulttype::Vec3d ;

class  SOFA_SOFAIMPLICITFIELD_API twist : public ScalarField
{
public:
    SOFA_CLASS(Twist, ScalarField);

public:
    twist() ;
    ~twist() override { }

    /// Inherited from BaseObject
    ///void init() override ;
    ///void reinit() override ;

    /// Inherited from ScalarField.
    ///double getValue(Vec3d& Pos, int &domain) override ;
    

    ///using ScalarField::getValue ;
    

    Data<bool> T_inside; ///< If true the field is oriented inside (resp. outside) the torus. (default = false)
    Data<double> T_sradius ; ///< The small radius of the torus. (default = 1)
    Data<double> T_lradius ; ///< The large radius of the torus. (default = 1)
     Data<Vec3d> T_center; ///< Position of the cube Surface. (default=0 0 0)

protected:
    double small_radius, large_radius;
    bool t_inside;
    Vec3d t_center;
};

} /// _torusfield_

using _torusfield_::TorusField ;

} /// geometry

} /// component

} /// sofa

#endif


