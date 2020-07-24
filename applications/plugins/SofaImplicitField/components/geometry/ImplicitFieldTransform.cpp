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

#include "ImplicitFieldTransform.h"

namespace sofa::core::objectmodel::BaseObject
{
///using sofa::defaulttype::RGBAColor;
///using sofa::defaulttype::Ray;
using sofa::defaulttype::Vec3d ;

/// Register in the Factory
static int ImplicitFieldTransformClass = core::RegisterObject("registering of ImplicitFieldTransform class") .add<ImplicitFieldTransform>();

//ImplicitFieldTransform::ImplicitFieldTransform() :
//    l_field(initLink("field", "The scalar field to render"))

Vec3d ImplicitFieldTransform::opTwist(Vec3d & P){
    double k = 3.0;
    double c = cos(k*P[1]);
    double s = sin(k*P[1]);
    Vec3d z ;
    z[0] = P[0]*c -P[2]*s;
    z[1] = P[0]*s +P[2]*c;
    z[2] = P[1];

    return z ;

}

Vec3d ImplicitFieldTransform::opcheapBend(Vec3d & P){
    double k = 3.0 ;
    double c = cos(k*P[1]);
    double s = sin(k*P[1]);
    Vec3d w ;
    w[0] = P[0]*c -P[1]*s;
    w[1] = P[0]*s +P[1]*c;
    w[2] = P[2];

    return w;
}


} /// namespace sofa::core::objectmodel::BaseObject


