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

#include <SofaImplicitField/config.h>
#include <sofa/core/objectmodel/BaseObject.h>
#include <SofaImplicitField/components/geometry/ScalarField.h>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/defaulttype/VecTypes.h>

namespace sofaimplicitfield
{

using sofa::core::objectmodel::BaseLink;
using sofa::core::objectmodel::SingleLink;
using sofa::component::geometry::ScalarField;
using sofa::defaulttype::Vec3d;
using sofa::defaulttype::Vec3Types;
using sofa::defaulttype::Vec4d;
//using sofa::defaulttype::Vec4Types;

class SOFA_SOFAIMPLICITFIELD_API DisplacementField : public sofaimplicitfield::ScalarField
{
public:
   SOFA_CLASS(DisplacementField, sofa::core::objectmodel::BaseObject);

   SingleLink<DisplacementField, ScalarField, BaseLink::FLAG_STRONGLINK> l_field;
   SingleLink<DisplacementField, sofa::core::topology::BaseMeshTopology, BaseLink::FLAG_STOREPATH> l_topology;
   SingleLink<DisplacementField, sofa::core::behavior::MechanicalState<Vec3Types>, BaseLink::FLAG_STOREPATH> l_dofs;

   int getDomain(Vec3d& pos, int domain) override;
   double getValue(Vec3d& pos, int& domain) override;
   Vec3d getGradient(Vec3d& pos, int& domain) override;

   Vec4d getBarycentricCoordinates(const Vec3d& p, int& domain, sofa::helper::ReadAccessor<sofa::helper::vector<Vec3d>>& dof);
   bool checkPointInTetrahedronAndGetBarycentricCoordinates(const Vec3d& p, int& domain, sofa::helper::ReadAccessor<sofa::helper::vector<Vec3d>>& dof, Vec4d& barycentric_coefs);

protected:
   DisplacementField();
   ~DisplacementField() override {}
   
   double determinant4x4ForVec3And1(const Vec3d& v0, const Vec3d& v1, const Vec3d& v2, const Vec3d& v3); // TODO: replace this by a call to BarycentricMapper?
   
   Vec4d getBarycentricCoordinates(const Vec3d& p, const Vec3d& v0, const Vec3d& v1, const Vec3d& v2, const Vec3d& v3);

   bool checkPointInTetrahedronAndGetBarycentricCoordinates(const Vec3d& p, const Vec3d& v0, const Vec3d& v1, const Vec3d& v2, const Vec3d& v3, Vec4d& barycentric_coefs);
};

}
