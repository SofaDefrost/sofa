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
int ImplicitFieldTransformClass = sofa::core::RegisterObject("registering of ImplicitFieldTransform class").add<DisplacementField>();

DisplacementField::DisplacementField() :
    l_field(initLink("field", "The scalar field to displace")),
    l_topology(initLink("topology", "The mesh topology to use as interpolation field")),
    l_dofs(initLink("dofs", "The nodal values to interpolate."))
{
/* TODO
    // Read DOFs current positions.
    auto dof = getReadAccessor(*l_dofs->read(sofa::core::VecCoordId::position()));
    // Read DOFs rest positions:
    auto dof_rest = getReadAccessor(*l_dofs->read(sofa::core::VecCoordId::restPosition()));
TODO */
}

/* overwritten */

double DisplacementField::getValue(Vec3d& pos, int& domain)
{
    SOFA_UNUSED(domain);
    // Initialise containers.
    bool found;
    Vec4d barycentric_coefs {0.0, 0.0, 0.0, 0.0};
    // Iterate over each tetrahedron:
    for (int t=0; t<l_topology->getNbTetrahedra(); t++)
    {
        // Test belonging and compute barycentric coefficients.
        found = checkPointInTetrahedronAndGetBarycentricCoordinates(pos, t, barycentric_coefs);
        if (found)
        {
            // Read DOFs rest positions:
            auto dof_rest = getReadAccessor(*l_dofs->read(sofa::core::VecCoordId::restPosition()));
            // Get target tetrahedron.
            auto tetra = l_topology->getTetrahedron(t);
            // Compute the underformed coordinate of 'pos':
            Vec3d pos_undeformed = barycentric_coefs[0] * dof_rest[tetra[0]] + 
                barycentric_coefs[1] * dof_rest[tetra[1]] + 
                barycentric_coefs[2] * dof_rest[tetra[2]] + 
                barycentric_coefs[3] * dof_rest[tetra[3]];
            return l_field->getValue(pos_undeformed, domain);
        }
    }
    // The point 'pos' doesn't belong to any tetrahedra. It isn't subject to a deformation.
    // msg_warning() << "The point " << pos << " does not belong to any tetrahedron of the displacement field!";
    return -1; // l_field->getValue(pos, domain);
}

Vec3d DisplacementField::getGradient(Vec3d& pos, int& domain)
{
    // Initialise containers.
    bool found;
    Vec3d gradient {0.0, 0.0, 0.0};
    Vec4d barycentric_coefs {0.0, 0.0, 0.0, 0.0};
    // Iterate over each tetrahedron:
    for ( int t=0; t<l_topology->getNbTetrahedra(); t++)
    {
        // Test belonging and compute barycentric coefficients.
        found = checkPointInTetrahedronAndGetBarycentricCoordinates(pos, t, barycentric_coefs);
        if (found)
        {
            // Evaluate point.
            double v = getValue(pos, t, domain);
            // Evaluate displaced point:
            double epsilon = d_epsilon.getValue();
            pos[0] += epsilon;
            gradient[0] = getValue(pos, t, domain);
            pos[0] -= epsilon;
            pos[1] += epsilon;
            gradient[1] = getValue(pos, t, domain);
            pos[1] -= epsilon;
            pos[2] += epsilon;
            gradient[2] = getValue(pos, t, domain);
            pos[2] -= epsilon;
            // Finite difference.
            gradient[0] = (gradient[0]-v)/epsilon;
            gradient[1] = (gradient[1]-v)/epsilon;
            gradient[2] = (gradient[2]-v)/epsilon;
            return gradient;
        }
    }
    // The point 'pos' doesn't belong to any tetrahedra. It isn't subject to a deformation.
    // msg_warning() << "The point " << pos << " does not belong to any tetrahedron of the displacement field!";
    return l_field->getGradient(pos, domain);
}

/* public */

int DisplacementField::getDeformationId(Vec3d& pos)
{
    // Initialize containers.
    bool found;
    Vec4d barycentric_coefs;
    // Iterate over each tetrahedron:
    for ( int t=0; t<l_topology->getNbTetrahedra(); t++)
    {
        // Test belonging and compute barycentric coefficients.
        found = checkPointInTetrahedronAndGetBarycentricCoordinates(pos, t, barycentric_coefs);
        if (found)
        {
            return t;
        }
    }
    return -1;
}


double DisplacementField::getValue(Vec3d& pos, int& tetrahedron_id, int& domain)
{
    SOFA_UNUSED(domain);
    // Read DOFs rest positions:
    auto dof_rest = getReadAccessor(*l_dofs->read(sofa::core::VecCoordId::restPosition()));
    // Get target tetrahedron.
    auto tetra = l_topology->getTetrahedron(tetrahedron_id);
    
/* Choose:
    // Initialise containers.
    bool found;
    Vec4d barycentric_coefs {0.0, 0.0, 0.0, 0.0};
    // Test belonging and compute barycentric coefficients.
    found = checkPointInTetrahedronAndGetBarycentricCoordinates(pos, tetrahedron_id, barycentric_coefs);
    if (!found)
    {
        // msg_warning() << "Attempting to compute the value at point (" << pos << ")  under the deformation of the "<< tetrahedron_id << "'s tetrahedron while it is outside said tetrahedron!";
        // TODO: Produce an error?
    }
    // Compute the barycentric coorinates of 'pos' in the tetrahedron:
Or: */
    Vec4d barycentric_coefs = getBarycentricCoordinates(pos, tetrahedron_id);

    // Compute the underformed coordinate of 'pos':
    Vec3d pos_undeformed = barycentric_coefs[0] * dof_rest[tetra[0]] + 
        barycentric_coefs[1] * dof_rest[tetra[1]] + 
        barycentric_coefs[2] * dof_rest[tetra[2]] + 
        barycentric_coefs[3] * dof_rest[tetra[3]];
    return l_field->getValue(pos_undeformed, domain);
}

Vec3d DisplacementField::getGradient(Vec3d& pos, int& tetrahedron_id, int& domain)
{
    // Initialise containers.
    Vec3d gradient {0.0, 0.0, 0.0};
    // Evaluate point.
    double v = getValue(pos, tetrahedron_id, domain);
    // Evaluate displaced point:
    double epsilon = d_epsilon.getValue();
    pos[0] += epsilon;
    gradient[0] = getValue(pos, tetrahedron_id, domain);
    pos[0] -= epsilon;
    pos[1] += epsilon;
    gradient[1] = getValue(pos, tetrahedron_id, domain);
    pos[1] -= epsilon;
    pos[2] += epsilon;
    gradient[2] = getValue(pos, tetrahedron_id, domain);
    pos[2] -= epsilon;
    // Finite difference.
    gradient[0] = (gradient[0]-v)/epsilon;
    gradient[1] = (gradient[1]-v)/epsilon;
    gradient[2] = (gradient[2]-v)/epsilon;
    return gradient;
}

Vec4d DisplacementField::getBarycentricCoordinates(const Vec3d& p, const int& tetrahedron_id)
{
    // Read DOFs current positions.
    auto dof = getReadAccessor(*l_dofs->read(sofa::core::VecCoordId::position()));
    // Get target tetrahedron.
    auto tetra = l_topology->getTetrahedron(tetrahedron_id);
    return getBarycentricCoordinates(p, dof[tetra[0]], dof[tetra[1]], dof[tetra[2]], dof[tetra[3]]);
}

bool DisplacementField::checkPointInTetrahedronAndGetBarycentricCoordinates(const Vec3d& p, const int& tetrahedron_id, Vec4d& barycentric_coefs)
{
    // Read DOFs current positions.
    auto dof = getReadAccessor(*l_dofs->read(sofa::core::VecCoordId::position()));
    // Get target tetrahedron.
    auto tetra = l_topology->getTetrahedron(tetrahedron_id);
    return checkPointInTetrahedronAndGetBarycentricCoordinates(p, dof[tetra[0]], dof[tetra[1]], dof[tetra[2]], dof[tetra[3]], barycentric_coefs);
}

/* protected */

/**
 * TODO: replace this by a call to BarycentricMapper?
 *
 * Compute the derterminant of the matrix:
 * ---                    ---
 * | v0[0], v0[1], v0[2], 1 |
 * | v1[0], v1[1], v1[2], 1 |
 * | v2[0], v2[1], v2[2], 1 |
 * | v3[0], v3[1], v3[2], 1 |
 * ---                    ---
 *
**/
double DisplacementField::determinant4x4ForVec3And1(const Vec3d& v0, const Vec3d& v1, const Vec3d& v2, const Vec3d& v3)
{
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

Vec4d DisplacementField::getBarycentricCoordinates(const Vec3d& p, const Vec3d& v0, const Vec3d& v1, const Vec3d& v2, const Vec3d& v3)
{
    Vec4d barycentric_coefs {0.0, 0.0, 0.0, 0.0};
    // Compute the tetrahedron's determinant.
    double d0 = determinant4x4ForVec3And1(v0, v1, v2, v3);
    // Compute the determinants of each sub-tetrahedra.
    double d1 = determinant4x4ForVec3And1(p, v1, v2, v3);
    double d2 = determinant4x4ForVec3And1(v0, p, v2, v3);
    double d3 = determinant4x4ForVec3And1(v0, v1, p, v3);
    double d4 = determinant4x4ForVec3And1(v0, v1, v2, p);
    // Compute the barycentric coeffcients.
    barycentric_coefs[0] = d1/d0;
    barycentric_coefs[1] = d2/d0;
    barycentric_coefs[2] = d3/d0;
    barycentric_coefs[3] = d4/d0;
    return barycentric_coefs;
}

bool DisplacementField::checkPointInTetrahedronAndGetBarycentricCoordinates(const Vec3d& p, const Vec3d& v0, const Vec3d& v1, const Vec3d& v2, const Vec3d& v3, Vec4d& barycentric_coefs)
{
    // Compute the tetrahedron's determinant.
    double d0 = determinant4x4ForVec3And1(v0, v1, v2, v3);
    // Compute the determinants of each sub-tetrahedra.
    double d1 = determinant4x4ForVec3And1(p, v1, v2, v3);
    double d2 = determinant4x4ForVec3And1(v0, p, v2, v3);
    double d3 = determinant4x4ForVec3And1(v0, v1, p, v3);
    double d4 = determinant4x4ForVec3And1(v0, v1, v2, p);
    if ((d0<0 && d1<=0 && d2<=0 && d3<=0 && d4<=0) || (d0>0 && d1>=0 && d2>=0 && d3>=0 && d4>=0))
    {
        // Compute the barycentric coeffcients.
        barycentric_coefs[0] = d1/d0;
        barycentric_coefs[1] = d2/d0;
        barycentric_coefs[2] = d3/d0;
        barycentric_coefs[3] = d4/d0;
        return true;
    }
    return false;
}

} /// sofaimplicitfield


