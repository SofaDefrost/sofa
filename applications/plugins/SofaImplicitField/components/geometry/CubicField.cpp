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

#include "CubicField.h"
#include<algorithm>

namespace sofa
{

namespace component
{

namespace geometry
{

namespace _cubicfield_
{
using sofa::defaulttype::Vec3d ;

CubicField::CubicField()
    : c_inside(initData(&c_inside, false, "inside", "If true the field is oriented inside (resp. outside) the cube. (default = false)"))
    ,c_edge(initData(&c_edge, 1.0, "edge", "side of cube emitting the field. (default = 1)"))
    , c_center(initData(&c_center, Vec3d(0.0,0.0,0.0), "center", "Position of the cubeSurface. (default=0 0 0)" ))
   
{init();
}

void CubicField::init()
{
    m_inside = c_inside.getValue();
    m_edge= c_edge.getValue();
    m_center = c_center.getValue();
    
}

void CubicField::reinit()
{
    init();
}

Vec3d CubicField::opTwist(Vec3d & P){
    double k = 3.0;
    double c = cos(k*P[1]);
    double s = sin(k*P[1]);
    //Matrix3d mat2 ;
//    defaulttype::Mat<2,2,double> mat2;
//    mat2[0][0] = c ;
//    mat2[0][1] = -s ;
//    mat2[1][0] = s;
//    mat2[1][1] = c ;
    Vec3d z ;
    z[0] = P[0]*c -P[2]*s;
    z[1] = P[0]*s +P[2]*c;
    z[2] = P[1];

    return z ;

}
Vec3d CubicField::opcheapBend(Vec3d & P){
    double k = 3.0 ;
    double c = cos(k*P[1]);
    double s = sin(k*P[1]);
//    defaulttype::Mat<2,2,double> mat2;
//    mat2[0][0] = c ;
//    mat2[0][1] = -s ;
//    mat2[1][0] = s;
//    mat2[1][1] = c ;
    Vec3d w ;
    w[0] = P[0]*c -P[1]*s;
    w[1] = P[0]*s +P[1]*c;
    w[2] = P[2];

    return w;
}

double CubicField::getValue(Vec3d& w, int& domain)
{
    SOFA_UNUSED(domain) ;
    Vec3d Pos = opTwist(w);
    Vec3d d ;
    Vec3d q ;
    /// Vec3d b ;
    double result ;
    d[0] = std::abs(Pos[0]) - m_edge;
    d[1] = std::abs(Pos[1]) - m_edge;
    d[2] = std::abs(Pos[2]) - m_edge;
    q[0] = std::max(d[0],0.0);
    q[1] = std::max(d[1],0.0);
    q[2] = std::max(d[2],0.0);
    double len_q = sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2]);

        result = std::max(len_q,0.0) + std::min(std::max(d[0],std::max(d[1],d[2])),0.0);

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

