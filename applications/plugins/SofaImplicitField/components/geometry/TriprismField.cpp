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
#include <sofa/core/ObjectFactory.h>
using sofa::core::RegisterObject ;

#include "TriprismField.h"
#include<algorithm>

namespace sofa
{

namespace component
{

namespace geometry
{

namespace _triprismfield_
{
using sofa::defaulttype::Vec3d ;
using sofa::defaulttype::Vec2d ;

TriprismField::TriprismField()
    : Tr_center(initData(&Tr_center, Vec3d(0.0,0.0,0.0), "center" , "height of the prism. default = 1 1"))
    , Tr_height(initData(&Tr_height, Vec2d(1.0,1.0), "height", "Position of the TriprismSurface. (default=0 0 0)" ))
   
{init();
}

void TriprismField::init()
{
    tr_center = Tr_center.getValue();
    h = Tr_height.getValue();
    
}

void TriprismField::reinit()
{
    init();

}

double TriprismField::getValue(Vec3d& Pos, int& domain)
{
    SOFA_UNUSED(domain) ;
    Vec3d q;
    q[0] = std::abs(Pos[0]);
    q[1] = std::abs(Pos[1]);
    q[2] = std::abs(Pos[2]);
    return std::max(q[0]-h[1],std::max(q[0]*0.866025+Pos[1]*0.5,-Pos[1])-h[0]*0.5);
}



// Register in the Factory
static int TorusFieldClass = core::RegisterObject("A Triprism implicit field.")
        .add< TriprismField >()
        ;

} /// _TriprismField_
} /// implicit
} /// component
} /// sofa



