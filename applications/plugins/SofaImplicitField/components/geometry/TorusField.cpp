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
#include <sofa/defaulttype/Mat.h>

using sofa::core::RegisterObject ;

#include "TorusField.h"
#include <Eigen/Sparse>
#include<algorithm>

namespace sofa
{

namespace component
{

namespace geometry
{

namespace _torusfield_
{
using sofa::defaulttype::Vec3d ;
using sofa::defaulttype::Vec2d ;
using namespace Eigen;


TorusField::TorusField()
    : T_inside(initData(&T_inside, false, "inside", "If true the field is oriented inside (resp. outside) the torus. (default = false)"))
    ,T_sradius(initData(&T_sradius, 0.5, "sradius", "the small radius of the torus. (default = 0.5)"))
    ,T_lradius(initData(&T_lradius, 1.0, "lradius", "the large radius of the torus. (default = 1)"))
    , T_center(initData(&T_center, Vec3d(0.0,0.0,0.0), "center", "Position of the TorusSurface. (default=0 0 0)" ))
   
{init();
}

void TorusField::init()
{
    t_inside = T_inside.getValue();
   small_radius= T_sradius.getValue();
  large_radius= T_lradius.getValue();
    t_center = T_center.getValue();
    
}

void TorusField::reinit()
{
    init();

}






Vec3d TorusField::opTwist(Vec3d & P){
    double k = 200.0 ;
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
double TorusField::getValue(Vec3d & Pos, int & domain)
{   ///msg_warning()<< "valzur avant twist" << w ;
    SOFA_UNUSED(domain);
     ///Vec3d Pos = opcheapBend(w);
    ///msg_warning()<< "valzur aprés twist" << Pos;
    double q = sqrt(Pos[0]*Pos[0] + Pos[1]*Pos[1]) - large_radius;
    double len_S = sqrt(q*q + Pos[2]*Pos[2]);
    double resultat = len_S- small_radius;
    return resultat;

}
double TorusField::opDisplace(Vec3d & P){
    double dist1 = getValue(P);
    double dist2 = sin(2.*P[0])*sin(2.*P[1])*sin(2.*P[2]);
    return dist1+dist2;
}

Vec3d TorusField::opcheapBend(Vec3d & P){
    double k = 1000.0 ;
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


// Register in the Factory
static int TorusFieldClass = core::RegisterObject("A Torus implicit field.")
        .add< TorusField >()
        ;

} /// _TorusField_
} /// implicit
} /// component
} /// sofa


