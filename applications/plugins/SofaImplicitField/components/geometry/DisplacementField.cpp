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
    
    // Read DOFs current positions.
    auto dof = getReadAccessor(*l_dofs->read(sofa::core::VecCoordId::position()));
    // Read DOFs rest positions:
    auto dof_rest = getReadAccessor(*l_dofs->read(sofa::core::VecCoordId::restPosition()));

    // Iterate over each tetrahedron:
    for (int t=0; l_topology->getNbTetrahedra(); t++)//(auto tetra : l_topology->getTetrahedra())
    {
        // Get current positions of the tetrahedron vertices.
        auto tetra = l_topology->getTetrahedron(t);
        Vec3d r0 = dof[tetra[0]];
        Vec3d r1 = dof[tetra[1]];
        Vec3d r2 = dof[tetra[2]];
        Vec3d r3 = dof[tetra[3]];

        /* determinant() is undefined for 4x4 matrices?
        Vec4d r0_1 = {r0[0], r0[1], r0[2], 1};
        etc...
        sofa::defaulttype::Mat4x4d T0 = {r0_1, r1_1, r2_1, r3_1};
        etc...
        auto d0 = determinant(T0);
        etc... */

        // Compute the determinants of each sub-tetrahedra.
        double d0 = determinant4x4ForVec3And1(r0, r1, r2, r3); // DisplacementField:: <- ?
        double d1 = determinant4x4ForVec3And1(pos, r1, r2, r3);
        double d2 = determinant4x4ForVec3And1(r0, pos, r2, r3);
        double d3 = determinant4x4ForVec3And1(r0, r1, pos, r3);
        double d4 = determinant4x4ForVec3And1(r0, r1, r2, pos);

        // Test if 'pos' belongs to the current tetrahedron.
        if (d0 == 0)
        {
            msg_warning() << "The tetrahedron is degenerate (flat)!";
            continue;
        }
        if ((d0<0 && d1<=0 && d2<=0 && d3<=0 && d4<=0) || (d0>0 && d1>=0 && d2>=0 && d3>=0 && d4>=0))
        {
            // Draws the deformed state.
            //drawtools->drawTetrahedron(r0, r1, r2, r3, RGBAColor::red());

            // Compute the barycentric coefficients of 'pos' in the deformed tetrahedron:
            double coef0 = d1/d0;
            double coef1 = d2/d0;
            double coef2 = d3/d0;
            double coef3 = d4/d0;

            // Compute the underformed coordinate of 'pos':
            Vec3d p0 = coef0 * dof_rest[tetra[0]] + coef1 * dof_rest[tetra[1]] + coef2 * dof_rest[tetra[2]] + coef3 * dof_rest[tetra[3]];
            
            msg_warning() << "The tetrahedron containing the point (" << pos << ") has been found!";
            msg_warning() << "d0: " << d0;
            msg_warning() << "d1: " << d1 << ", coef0: " << coef0;
            msg_warning() << "d2: " << d2 << ", coef1: " << coef1;
            msg_warning() << "d3: " << d3 << ", coef2: " << coef2;
            msg_warning() << "d4: " << d4 << ", coef3: " << coef3;
            msg_warning() << "new position: " << p0;
            msg_warning() << "l_field->getValue( p0 ): " << l_field->getValue( p0 ) << "\n";

            return l_field->getValue( p0 );
        }
    }
    // The point 'pos' doesn't belong to any tetrahedra. It isn't subject to a deformation.
    msg_warning() << "The point " << pos << " does not belong to any tetrahedron of the displacement field!";
    return l_field->getValue(pos);
}

// TODO: replace this by a call to BarycentricMapper?
double DisplacementField::determinant4x4ForVec3And1(Vec3d& v0, Vec3d& v1, Vec3d& v2, Vec3d& v3)
{
    /*
    Compute the derterminant of the matrix:
    ---                    ---
    | v0[0], v0[1], v0[2], 1 |
    | v1[0], v1[1], v1[2], 1 |
    | v2[0], v2[1], v2[2], 1 |
    | v3[0], v3[1], v3[2], 1 |
    ---                    ---
    */
    double det = v1[2]*v2[1]*v3[0] - v0[2]*v2[1]*v3[0] -
        v1[1]*v2[2]*v3[0] + v0[1]*v2[2]*v3[0] +
        v0[2]*v1[1]*v3[0] - v0[1]*v1[2]*v3[0] -
        v1[2]*v2[0]*v3[1] + v0[2]*v2[0]*v3[1] +
        v1[0]*v2[2]*v3[1] - v0[0]*v2[2]*v3[1] -
        v0[2]*v1[0]*v3[1] + v0[0]*v1[2]*v3[1] +
        v1[1]*v2[0]*v3[2] - v0[1]*v2[0]*v3[2] -
        v1[0]*v2[1]*v3[2] + v0[0]*v2[1]*v3[2] +
        v0[1]*v1[0]*v3[2] - v0[0]*v1[1]*v3[2] -
        v0[2]*v1[1]*v2[0] + v0[1]*v1[2]*v2[0] +
        v0[2]*v1[0]*v2[1] - v0[0]*v1[2]*v2[1] -
        v0[1]*v1[0]*v2[2] + v0[0]*v1[1]*v2[2];
    return det;
}

} /// sofaimplicitfield


