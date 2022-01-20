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
#include <sofa/core/objectmodel/Data.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/type/RGBAColor.h>
#include <sofa/gl/GLSLShader.h>

namespace sofaimplicitfield
{

using sofa::helper::getReadAccessor;
using sofa::type::RGBAColor;
using sofa::type::Mat3x3d;

/// Register in the Factory
int ImplicitFieldTransformClass = sofa::core::RegisterObject("registering of ImplicitFieldTransform class").add<DisplacementField>();

<<<<<<< HEAD
class DisplacementField::InternalData
{
public:
    GLuint vertexBufferObject;
    GLuint displacementBufferObject;
    GLuint barycentriCoordinatesBufferObject;

    std::vector<sofa::type::Vec3f> vertices;
    std::vector<sofa::type::Vec3f> displacements;
    std::vector<unsigned int> indices;
    bool isInited = false;
    GLuint mvpMatrixID ;

    void init()
    {
        std::cout << "DisplacementField init internal data" << std::endl;
        isInited = true;

        /// Create a vertex buffer object so we can draw it
        glGenBuffers(1, &vertexBufferObject);

        /// Create a vertex buffer object so we can draw it
        glGenBuffers(1, &displacementBufferObject);

        /// Create a vertex buffer object so we can draw it
        glGenBuffers(1, &barycentriCoordinatesBufferObject);
    }
};

DisplacementField::DisplacementField() :
=======
DisplacementField::DisplacementField():
>>>>>>> defrost/stage-nicola-zotto
    l_field(initLink("field", "The scalar field to displace")),
    l_topology(initLink("topology", "The mesh topology to use as interpolation field")),
    l_dofs(initLink("dofs", "The nodal values to interpolate.")),
    l_shader(initLink("shader", "The shader to use for the rendering."))
{
    data.reset(new InternalData());
}

/* overwritten */

double DisplacementField::getValue(Vec3d& pos, int& domain)
{
    // Initialise containers.
    bool found;
    Vec4d barycentric_coefs {0.0, 0.0, 0.0, 0.0};
    // Initialize accessors.
    // Read DOFs rest positions.
    auto dof_rest = getReadAccessor(*l_dofs->read(sofa::core::VecCoordId::restPosition()));
    // Read DOFs current positions.
    auto dof = getReadAccessor(*l_dofs->read(sofa::core::VecCoordId::position()));
    // If a domain was specified, check validity.
    if (domain!=-1)
    {
        // Test belonging and compute barycentric coefficients.
        found = checkPointInTetrahedronAndGetBarycentricCoordinates(pos, domain, dof, barycentric_coefs);
        if (found)
        {
            // Get target tetrahedron.
            auto tetra = l_topology->getTetrahedron(domain);
            // Compute the underformed coordinate of 'pos':
            Vec3d pos_undeformed = barycentric_coefs[0] * dof_rest[tetra[0]] + 
                barycentric_coefs[1] * dof_rest[tetra[1]] + 
                barycentric_coefs[2] * dof_rest[tetra[2]] + 
                barycentric_coefs[3] * dof_rest[tetra[3]];
            return l_field->getValue(pos_undeformed, domain);
        }
        // msg_warning() << "The point " << pos << " was not found in the passed domain: << domain << "!";
    }
    // Iterate over each tetrahedron:
    for (int t=0; t<l_topology->getNbTetrahedra(); t++)
    {
        // Test belonging and compute barycentric coefficients.
        found = checkPointInTetrahedronAndGetBarycentricCoordinates(pos, t, dof, barycentric_coefs);
        if (found)
        {
            // Save found domain.
            domain = t;
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
    // msg_warning() << "The point " << pos << " does not belong to any tetrahedron of the displacement field!";
    // Forget/Reset domain.
    domain = -1;
    // The point 'pos' doesn't belong to any tetrahedra. It isn't subject to a deformation.
    return l_field->getValue(pos, domain);
}

Vec3d DisplacementField::getGradient(Vec3d& pos, int& domain)
{
    // Initialise containers.
    bool found;
    Vec3d gradient {0.0, 0.0, 0.0};
    Vec4d barycentric_coefs {0.0, 0.0, 0.0, 0.0};
    // Initialize accessors.
    // Read DOFs current positions.
    auto dof = getReadAccessor(*l_dofs->read(sofa::core::VecCoordId::position()));
    // If a domain was specified, check validity.
    if (domain!=-1)
    {
        // Test belonging and compute barycentric coefficients.
        found = checkPointInTetrahedronAndGetBarycentricCoordinates(pos, domain, dof, barycentric_coefs);
        if (found)
        {
            // Evaluate point.
            double v = getValue(pos, domain);
            // Evaluate displaced point:
            double epsilon = d_epsilon.getValue(); // l_field->d_epsilon.getValue();
            pos[0] += epsilon;
            gradient[0] = getValue(pos, domain);
            pos[0] -= epsilon;
            pos[1] += epsilon;
            gradient[1] = getValue(pos, domain);
            pos[1] -= epsilon;
            pos[2] += epsilon;
            gradient[2] = getValue(pos, domain);
            pos[2] -= epsilon;
            // Finite difference.
            gradient[0] = (gradient[0]-v)/epsilon;
            gradient[1] = (gradient[1]-v)/epsilon;
            gradient[2] = (gradient[2]-v)/epsilon;
            return gradient;
        }
        // msg_warning() << "The point " << pos << " was not found in the passed domain: << domain << "!";
    }
    // Iterate over each tetrahedron:
    for (int t=0; t<l_topology->getNbTetrahedra(); t++)
    {
        // Test belonging and compute barycentric coefficients.
        found = checkPointInTetrahedronAndGetBarycentricCoordinates(pos, t, dof, barycentric_coefs);
        if (found)
        {
            // Save found domain.
            domain = t;
            // Evaluate point.
            double v = getValue(pos, t);
            // Evaluate displaced point:
            double epsilon = d_epsilon.getValue(); //l_field->d_epsilon.getValue();
            pos[0] += epsilon;
            gradient[0] = getValue(pos, t);
            pos[0] -= epsilon;
            pos[1] += epsilon;
            gradient[1] = getValue(pos, t);
            pos[1] -= epsilon;
            pos[2] += epsilon;
            gradient[2] = getValue(pos, t);
            pos[2] -= epsilon;
            // Finite difference.
            gradient[0] = (gradient[0]-v)/epsilon;
            gradient[1] = (gradient[1]-v)/epsilon;
            gradient[2] = (gradient[2]-v)/epsilon;
            return gradient;
        }
    }
    // msg_warning() << "The point " << pos << " does not belong to any tetrahedron of the displacement field!";
    // Forget/Reset domain.
    domain = -1;
    // The point 'pos' doesn't belong to any tetrahedra. It isn't subject to a deformation.
    return l_field->getGradient(pos, domain);
}

int DisplacementField::getDomain(Vec3d& pos, int domain)
{
    SOFA_UNUSED(domain);
    // Initialize containers.
    bool found;
    Vec4d barycentric_coefs;
    // Initialize accessors.
    // Read DOFs current positions.
    auto dof = getReadAccessor(*l_dofs->read(sofa::core::VecCoordId::position()));
    // Iterate over each tetrahedron:
    for (int t=0; t<l_topology->getNbTetrahedra(); t++)
    {
        // Test belonging and compute barycentric coefficients.
        found = checkPointInTetrahedronAndGetBarycentricCoordinates(pos, t, dof, barycentric_coefs);
        if (found)
        {
            return t;
        }
    }
    return -1;
}

<<<<<<< HEAD
    // Here are the tetrahedron's descriptions
    // l_topology->

    // Here are the moving position.
    // l_dofs->
=======
/* public */

Vec4d DisplacementField::getBarycentricCoordinates(const Vec3d& p, int& domain, sofa::helper::ReadAccessor<sofa::helper::vector<Vec3d>>& dof)
{
    // Get target tetrahedron.
    auto tetra = l_topology->getTetrahedron(domain);
    return getBarycentricCoordinates(p, dof[tetra[0]], dof[tetra[1]], dof[tetra[2]], dof[tetra[3]]);
}

bool DisplacementField::checkPointInTetrahedronAndGetBarycentricCoordinates(const Vec3d& p, int& domain, sofa::helper::ReadAccessor<sofa::helper::vector<Vec3d>>& dof, Vec4d& barycentric_coefs)
{
    // Get target tetrahedron.
    auto tetra = l_topology->getTetrahedron(domain);
    return checkPointInTetrahedronAndGetBarycentricCoordinates(p, dof[tetra[0]], dof[tetra[1]], dof[tetra[2]], dof[tetra[3]], barycentric_coefs);
}
>>>>>>> defrost/stage-nicola-zotto

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

<<<<<<< HEAD
/// Debug rendering of the Displacement Field.
void DisplacementField::draw(const sofa::core::visual::VisualParams* v)
{
    auto x = getReadAccessor(*l_dofs->read(sofa::core::VecCoordId::position()));
    auto x_0 = getReadAccessor(*l_dofs->read(sofa::core::VecCoordId::restPosition()));
    auto dt = v->drawTool();

    for(auto tetra : l_topology->getTetrahedra())
    {
        Vec3d r0 = x[tetra[0]];
        Vec3d r1 = x[tetra[1]];
        Vec3d r2 = x[tetra[2]];
        Vec3d r3 = x[tetra[3]];

        // Draws the initial state
        Mat3x3d T = {r0-r3, r1-r3, r2-r3};
        T.transpose();
        Mat3x3d Tinv = T.inverted();

        Vec3d coef0 = Tinv * (r0-r3);
        Vec3d coef1 = Tinv * (r1-r3);
        Vec3d coef2 = Tinv * (r2-r3);
    }

    if(!data.get())
        return;

    if(!data->isInited)
        data->init();

    data->vertices.clear();
    data->displacements.clear();
    for(sofa::Index i=0;i<x.size();++i)
    {
        auto& vertex = x[i];
        auto& vertex_at_rest_position = x_0[i];
        data->vertices.push_back(vertex);
        data->displacements.push_back( vertex_at_rest_position );
    }

    std::vector<unsigned int> indices;
    for(auto tetra : l_topology->getTetrahedra())
    {
        indices.push_back(tetra[0]);
        indices.push_back(tetra[1]);
        indices.push_back(tetra[2]);

        indices.push_back(tetra[1]);
        indices.push_back(tetra[0]);
        indices.push_back(tetra[3]);

        indices.push_back(tetra[0]);
        indices.push_back(tetra[2]);
        indices.push_back(tetra[3]);

        indices.push_back(tetra[2]);
        indices.push_back(tetra[1]);
        indices.push_back(tetra[3]);
    }

    ////// Compute sphere and depth
    double projMat[16];
    double modelMat[16];

    v->getProjectionMatrix(projMat);
    float fProjMat[16];
    for (unsigned int i = 0; i < 16; i++)
        fProjMat[i] = float(projMat[i]);

    v->getModelViewMatrix(modelMat);
    float fModelMat[16];
    for (unsigned int i = 0; i < 16; i++)
        fModelMat[i] = float(modelMat[i]);

    if(l_shader.get())
    {
        auto pm = l_shader->getUniform(l_shader->getCurrentIndex(), "projection_matrix");
        auto mm = l_shader->getUniform(l_shader->getCurrentIndex(), "object_matrix");
        auto uZNear = l_shader->getUniform(l_shader->getCurrentIndex(), "zNear");
        auto uZFar = l_shader->getUniform(l_shader->getCurrentIndex(), "zFar");

        l_shader->start();
        glUniform1f(uZNear, v->zNear());
        glUniform1f(uZFar, v->zFar());

        glUniformMatrix4fv(pm, 1, false, fProjMat);
        glUniformMatrix4fv(mm, 1, false, fModelMat);
    }else{
        msg_error() << "NO SHADER for rendering";
    }

    //glEnable(GL_DEPTH_TEST);
    //glDepthFunc(GL_LESS);

    //// Upload the geometry to the buffer.
    glBindBuffer(GL_ARRAY_BUFFER, data->vertexBufferObject);
    glBufferData(GL_ARRAY_BUFFER,                                   // target, size in byte of the buffer, buffer adress,
                 sizeof(sofa::type::Vec3f)*data->vertices.size(),
                 data->vertices.data(), GL_DYNAMIC_DRAW);

    //// Upload the displacement data into the buffer.
    glBindBuffer(GL_ARRAY_BUFFER, data->displacementBufferObject);
    glBufferData(GL_ARRAY_BUFFER,
                 sizeof(sofa::type::Vec3f)*data->displacements.size(),
                 data->displacements.data(), GL_DYNAMIC_DRAW);

    /// Position
    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, data->vertexBufferObject);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);

    /// Displacement
    glEnableVertexAttribArray(1);                                   // activates the generic vertex array
    glBindBuffer(GL_ARRAY_BUFFER, data->displacementBufferObject);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, 0);          // define an array of generic vertex attribute data

    /// Type, count, format, pointer
    glDrawElements(GL_TRIANGLES, 100, GL_UNSIGNED_INT, indices.data());
    glDisableVertexAttribArray(0);
    glDisableVertexAttribArray(1);
    if(l_shader.get())
        l_shader->stop();

}

} /// sofaimplicitfield
=======
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
>>>>>>> defrost/stage-nicola-zotto

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
