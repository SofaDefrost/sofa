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


#include <SofaImplicitField/config.h>
#include <sofa/core/objectmodel/BaseObject.h>
#include <SofaImplicitField/components/geometry/ScalarField.h>
///#include <SofaBaseVisual/BaseCamera.h>
///#include <sofa/core/DataEngine.h>



///using sofa::component::visualmodel::BaseCamera;
using sofa::component::geometry::ScalarField;
///using sofa::core::DataEngine;
using sofa::defaulttype::Vec3d ;
///using sofa::defaulttype::Vec2i;

class SOFA_SOFAIMPLICITFIELD_API ImplicitFieldTransform : public sofa::core::objectmodel::BaseObject
{
public:
   SOFA_CLASS(ImplicitFieldTransform, sofa::core::objectmodel::BaseObject);

   ///void init() override;
   ///void reinit() override;
   Vec3d opTwist(Vec3d & P);
   Vec3d  opcheapBend(Vec3d &P);
   ///void doUpdate() override;

   SingleLink<ImplicitFieldTransform, ScalarField, BaseLink::FLAG_STRONGLINK> l_field;

protected:
   ImplicitFieldTransform();
   ~ImplicitFieldTransform() override {}

  ///Vec3d transform(Vec3d (*operator)(Vec3d &P);
};

