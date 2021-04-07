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

#include "cubefield.h"

namespace sofa
{

namespace component
{

namespace geometry
{

namespace _Cubicfield_
{

CubicField::CubicField()
    : c_inside(initData(&d_inside, false, "inside", "If true the field is oriented inside (resp. outside) the cube. (default = false)"))
    , c_side(initData(&c_side, 1.0, "side", "side of the cube emitting the field. (default = 1)"))
   
{init();
}

void CubicField::init()
{
    m_inside = c_inside.getValue();
    b = c_side.getValue();
    
}

void CubicField::reinit()
{
    init();
}

double CubicField::getValue(Vec3d& Pos, int& domain)
{
    SOFA_UNUSED(domain) ;
    vec3 d ;
    d[0] = abs(pos[0]) - b[0];
    d[1] = abs(pos[1]) - b[1];
    d[2] = abs(pos[2]) - b[2];
    double len_d = sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
    if(m_inside)
        double result = max(len_d,0);
    else{
    	result = min(max(d[0],max(d[1],d[2])),0);
}
	 return result;
}



// Register in the Factory
static int CubicFieldClass = core::RegisterObject("A cubic implicit field.")
        .add< CubicField >()
        ;

} /// _CubicField_
} /// implicit
} /// component
} /// sofa

