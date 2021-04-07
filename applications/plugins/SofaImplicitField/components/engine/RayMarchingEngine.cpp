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

#include "RayMarchingEngine.h"

namespace sofa::component::engine
{
using sofa::defaulttype::RGBAColor;
using sofa::defaulttype::Ray;

/// Register in the Factory
static int RayMarchingEngineEngineClass = core::RegisterObject("A software implementation of the ray marching algorithm. "
                                                               "The method trace rays starting from the camera to locate a distance field. "
                                                               "The intersecting points are returned in a list of points")
        .add< RayMarchingEngine >();


RayMarchingEngine::RayMarchingEngine() :
    l_field(initLink("field", "The scalar field to render")),
    d_intersections(initData(&d_intersections, "intersections", "The intersection locations")),
    d_resolution(initData(&d_resolution, Vec2i{20,20}, "resolution", "The amount of samples per visual axis.")),
    d_renderingCube(initData(&d_renderingCube, {-2.0,2.0,-2.0,2.0,-2.0,2.0}, "renderingCube", "The cube from where the rendering is done"))
{
    addOutput(&d_intersections);
}

void RayMarchingEngine::init()
{
    setDirtyValue();
}

void RayMarchingEngine::reinit()
{
    setDirtyValue();
    update();
}

bool RayMarchingEngine::raymarch(const Ray& r, Vec3& out, ScalarField* field)
{
    double step=0.01;
    double cdist = 0.0;
    Vec3 cp = r.origin();
    const Vec3& dir = r.direction();
    for(unsigned int i=0;i<10 || cdist < 20; i++)
    {
        double d = field->getValue(cp);
        //std::cout << "Eval.." << cp << " => " << d << ", " << cdist << std::endl;
        if(fabs(d)<0.01)
        {
            out = cp;
            return true;
        }
        cdist += fabs(step);
        cp = cp + dir * fabs(step);
    }
    return false;
}

void RayMarchingEngine::doUpdate()
{
    std::cout << "Trying to updating the ray marching..." << std::endl;

    if(l_field.empty())
        return;

    ///// THEN tell everthing is (will be) up to date now
    /// @warning This must be done AFTER updating all inputs
    /// can be done before or after setting up the outputs
    cleanDirty();

    ScalarField* field = l_field.get();

    std::cout << "Updating the ray marching..." << std::endl;
    auto intersections = sofa::helper::getWriteAccessor(d_intersections);
    auto box = sofa::helper::getReadAccessor(d_renderingCube);
    auto res = sofa::helper::getReadAccessor(d_resolution);
    double rcw = res->x();
    double rch = res->y();
    double spacingw = 1.0 / rcw;
    double spacingh = 1.0 / rch;

    intersections.clear();
    Vec3 minb=box->minBBox();
    Vec3 maxb=box->maxBBox();
    Vec3 diffb=maxb-minb;
    std::vector<sofa::defaulttype::Vec2i> t {{1,2}, {1,2},
                                             {0,1}, {0,1},
                                             {0,2}, {0,2}};
    std::vector<Vec3> begin{minb, {maxb.x(), minb.y(), minb.z()},
                            minb, {minb.x(), minb.y(), maxb.z()},
                            minb, {minb.x(), maxb.y(), minb.z()}
                           };
    std::vector<double> d {1.0,-1.0, 1.0,-1.0, 1.0,-1.0};

    /// We iterate over the vectors...
    for(size_t k=0;k<d.size();k++)
    {
        auto& a = t[k];
        /// Compute the iteratoring parameter
        Vec3 nextline {0.0,0.0,0.0};
        Vec3 nextcolumn {0.0,0.0,0.0};
        nextcolumn[a[0]] = spacingw * diffb[a[0]];
        nextline[a[1]] = spacingh * diffb[a[1]];

        /// We generate the directional vector
        Vec3 dir1  {d[k],d[k],d[k]};
        dir1[a[0]] = 0.0;
        dir1[a[1]] = 0.0;

        Vec3 currL=begin[k];
        Vec3 out;
        for(unsigned int i = 0; i < rcw; i++)
        {
            Vec3 currP = currL;
            for(unsigned int j = 0; j < rch; j++)
            {
                Ray ray { currP, dir1 };
                if( raymarch(ray, out, field) )
                {
                    intersections.push_back(out);
                }
                currP += nextcolumn;
            }
            currL += nextline;
        }
    }
}

void RayMarchingEngine::draw(const sofa::core::visual::VisualParams* params)
{
    if(!params->displayFlags().getShowBehavior())
        return;
    auto dt = params->drawTool();
    double rcw = 40;
    double rch = 40;
    double spacingw = 1.0 / rcw;
    double spacingh = 1.0 / rch;

    auto box = sofa::helper::getReadAccessor(d_renderingCube);

    Vec3 minb=box->minBBox();
    Vec3 maxb=box->maxBBox();
    Vec3 diffb=maxb-minb;

    std::vector<sofa::defaulttype::Vec2i> t {{1,2}, {1,2},
                                             {0,1}, {0,1},
                                             {0,2}, {0,2}};
    std::vector<Vec3> begin{minb, {maxb.x(), minb.y(), minb.z()},
                            minb, {minb.x(), minb.y(), maxb.z()},
                            minb, {minb.x(), maxb.y(), minb.z()}
                           };
    std::vector<double> d {1.0,-1.0, 1.0,-1.0, 1.0,-1.0};

    /// We iterate over the vectors...
    for(size_t k=0;k<d.size();k++)
    {
        auto& a = t[k];

        /// Compute the iteratoring parameter
        Vec3 nextline {0.0,0.0,0.0};
        Vec3 nextcolumn {0.0,0.0,0.0};
        nextcolumn[a[0]] = spacingw * diffb[a[0]];
        nextline[a[1]] = spacingh * diffb[a[1]];

        /// We generate the directional vector
        Vec3 dir1  {d[k],d[k],d[k]};
        dir1[a[0]] = 0.0;
        dir1[a[1]] = 0.0;

        Vec3 currL = begin[k];
        Vec3 out;
        for(unsigned int i = 0; i < rcw; i++)
        {
            Vec3 currP = currL;
            for(unsigned int j = 0; j < rch; j++)
            {
                dt->drawSphere(currP, 0.02, RGBAColor::white());
                currP += nextcolumn;
                dt->drawLine(currP, currP+dir1*0.2, RGBAColor::blue());
            }
            currL += nextline;
        }
    }
}


} /// namespace sofa::component::engine

