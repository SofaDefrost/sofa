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

namespace _cubicfield_
{

using sofa::defaulttype::Vec3d ;

class  SOFA_SOFAIMPLICITFIELD_API CubicField  : public ScalarField
{
public:
    SOFA_CLASS(CubicField, ScalarField);

public:
    CubicField() ;
    ~CubicField() override { }

    /// Inherited from BaseObject
    void init() override ;
    void reinit() override ;

    /// Inherited from ScalarField.
    double getValue(Vec3d& Pos, int &domain) override ;
    Vec3d opTwist(Vec3d & P);
    Vec3d  opcheapBend(Vec3d &P);
    

    using ScalarField::getValue ;
    

    Data<bool> c_inside; ///< If true the field is oriented inside (resp. outside) the cube. (default = false)
    Data<double> c_edge ; ///< edge of the cube emitting the field. (default = 1)
     Data<Vec3d> c_center; ///< Position of the cube Surface. (default=0 0 0)

protected:
    double m_edge;
    bool m_inside;
    Vec3d m_center;
};

} /// _cubefield_

using _cubicfield_::CubicField ;

} /// geometry

} /// component

} /// sofa

#endif
