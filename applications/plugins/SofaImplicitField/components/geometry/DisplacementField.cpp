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

#include "DisplacementField.h"
#include <sofa/core/objectmodel/Data.h>

namespace sofaimplicitfield
{

using sofa::helper::getReadAccessor;

/// Register in the Factory
int ImplicitFieldTransformClass = sofa::core::RegisterObject("registering of ImplicitFieldTransform class") .add<DisplacementField>();

DisplacementField::DisplacementField() :
    l_field(initLink("field", "The scalar field to displace")),
    l_topology(initLink("topology", "The mesh topology to use as interpolation field")),
    l_dofs(initLink("dofs", "The nodal values to interpolate."))
{
}

double DisplacementField::getValue(Vec3d &pos, int &domain)
{
    SOFA_UNUSED(domain);

    /// Here are the tetrahedron's descriptions
    /// l_topology->

    /// Here are the moving position.
    /// l_dofs->

    auto x = getReadAccessor(*l_dofs->read(sofa::core::VecCoordId::position()));

    Vec3d p0 = x[0] - pos;
    return l_field->getValue( p0 );
}

} /// sofaimplicitfield


