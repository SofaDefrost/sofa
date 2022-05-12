/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, development version     *
*                (c) 2006-2017 INRIA, USTL, UJF, CNRS, MGH                    *
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
#define SOFA_COMPONENT_MAPPING_IMPLICITSURFACEMAPPING_CPP
#include <sofa/core/visual/VisualParams.h>
using sofa::core::visual::VisualParams ;

#include <sofa/core/ObjectFactory.h>
using sofa::core::RegisterObject ;

using sofa::type::Vec3dTypes ;
// using sofa::defaulttype::ExtVec3fTypes ;
using sofa::defaulttype::Vec4f ;

#include <sofa/helper/MarchingCubeUtility.h>

#include "SurfaceMeshGenerationFromScalarField.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
namespace sofa
{

namespace component
{

namespace geometry
{

namespace _surfacemeshgenerationfromscalarfield_
{

/// Register in the Factory
int SurfaceMeshGenerationFromScalarFieldClass =
        RegisterObject("Generate a surface mesh from an implicit function")
        .add< SurfaceMeshGenerationFromScalarField >() ;


SurfaceMeshGenerationFromScalarField::SurfaceMeshGenerationFromScalarField()
    : l_field(initLink("field", "The scalar field to generate a mesh from.")),
      mStep(initData(&mStep,0.1,"step","Step")),
      mIsoValue(initData(&mIsoValue,0.0,"isoValue","Iso Value")),
      mGridMin(initData(&mGridMin, Vec3d(-1,-1,-1),"min","Grid Min")),
      mGridMax(initData(&mGridMax, Vec3d(1,1,1),"max","Grid Max"))
      , d_out_points(initData(&d_out_points, "outputPoints", "position of the tiangles vertex"))
      , d_out_triangles(initData(&d_out_triangles, "outputTriangles", "list of triangles"))
      , d_visuCavity(initData(&d_visuCavity,true, "visuCavity","Allows to display the cavity surface"))
{
}

SurfaceMeshGenerationFromScalarField::~SurfaceMeshGenerationFromScalarField()
{
}

void SurfaceMeshGenerationFromScalarField::init()
{
    generateSurfaceMesh() ;
}

void SurfaceMeshGenerationFromScalarField::reinit()
{
}
static bool inited = false ;

void SurfaceMeshGenerationFromScalarField::draw(const VisualParams* vparams)
{


    if(!inited){
        //if( l_field.get() != 0 ){
            generateSurfaceMesh();
        //}
    }

    if(!d_visuCavity.getValue())
    {
        return;
    }

    vparams->drawTool()->drawBoundingBox(mGridMin.getValue(), mGridMax.getValue()) ;

    sofa::helper::ReadAccessor< Data<VecCoord> > x = d_out_points;
    sofa::helper::ReadAccessor< Data<SeqTriangles> > triangles = d_out_triangles;
    vparams->drawTool()->setLightingEnabled(true);

    for(const Triangle& triangle : triangles)
    {
        int a = triangle[0];
        int b = triangle[1];
        int c = triangle[2];
        Vec3d center = (x[a]+x[b]+x[c])*0.333333;
        Vec3d pa = (0.9*x[a]+0.1*center) ;
        Vec3d pb = (0.9*x[b]+0.1*center) ;
        Vec3d pc = (0.9*x[c]+0.1*center) ;

        Vec3d a1 = x[c]-x[b] ;
        Vec3d a2 = x[a]-x[b] ;
        
        vparams->drawTool()->drawTriangles({pa,pb,pc},
                                           {a1.cross(a2)},
                                           sofa::type::RGBAColor(0.0,0.0,1.0,1.0));
    }

    vparams->drawTool()->drawSpheres(d_out_points.getValue(), 0.01, sofa::type::RGBAColor(1.0,1.0,0.0,0.2)) ;
}

void SurfaceMeshGenerationFromScalarField::generateSurfaceMesh()
{
    if( l_field.get() == 0 )
        return ;

    if( !l_field.get()->isComponentStateValid() )
    {
        std::cout << "INVALID COMPONENT" << std::endl;
        return ;
    }
    inited=true;

    /// Copy the surface as a Sofa topology
    sofa::helper::WriteOnlyAccessor< Data<VecCoord> >      points = d_out_points;
    sofa::helper::WriteOnlyAccessor< Data<SeqTriangles> >  triangles = d_out_triangles;

    points.clear();
    triangles.clear();

    double isoval = mIsoValue.getValue();
    double mstep = mStep.getValue();
    double invStep = 1.0/mStep.getValue();

    Vec3d gridmin = mGridMin.getValue() ;
    Vec3d gridmax = mGridMax.getValue() ;

    int nx = floor((gridmax.x() - gridmin.x()) * invStep) + 1 ;
    int ny = floor((gridmax.y() - gridmin.y()) * invStep) + 1 ;
    int nz = floor((gridmax.z() - gridmin.z()) * invStep) + 1 ;

    double cx,cy,cz;
    int x,y,z,i,mk;
    const int *tri;

    (*planes.beginEdit()).resize(2*(nx)*(ny));
    P0 = (*planes.beginEdit()).begin()+0;
    P1 = (*planes.beginEdit()).begin()+nx*ny;

    const int dx = 1;
    const int dy = nx;

    z = 0;
    newPlane();

    i = 0 ;
    cz = gridmin.z()  ;
    for (int y = 0 ; y < ny ; ++y)
    {
        cy = gridmin.y() + mstep * y ;
        for (int x = 0 ; x < nx ; ++x, ++i)
        {
            cx = gridmin.x() + mstep * x ;

            Vec3d pos { cx, cy, cz }  ;
            double res =  l_field.get()->getValue(pos) ;
            (P1+i)->data = res ;
        }
    }


    for (z=1; z<=nz; ++z)
    {

        newPlane();

        i = 0 ;
        cz = gridmin.z() + mstep * z ;
        for (int y = 0 ; y < ny ; ++y)
        {
            cy = gridmin.y() + mstep * y ;
            for (int x = 0 ; x < nx ; ++x, ++i)
            {
                cx = gridmin.x() + mstep * x ;

                Vec3d pos { cx, cy, cz }  ;
                double res =  l_field.get()->getValue(pos) ;
                (P1+i)->data = res ;
            }
        }

        unsigned int i=0;
        int edgecube[12];
        const int edgepts[12] = {0,1,0,1,0,1,0,1,2,2,2,2};
        typename std::vector<CubeData>::iterator base = (*planes.beginEdit()).begin();
        int ip0 = P0-base;
        int ip1 = P1-base;
        edgecube[0]  = (ip0   -dy);
        edgecube[1]  = (ip0      );
        edgecube[2]  = (ip0      );
        edgecube[3]  = (ip0-dx   );
        edgecube[4]  = (ip1   -dy);
        edgecube[5]  = (ip1      );
        edgecube[6]  = (ip1      );
        edgecube[7]  = (ip1-dx   );
        edgecube[8]  = (ip1-dx-dy);
        edgecube[9]  = (ip1-dy   );
        edgecube[10] = (ip1      );
        edgecube[11] = (ip1-dx   );

        // First line is all zero
        {
            y=0;
            x=0;
            i+=nx;
        }
        for(y=1; y<ny; y++)
        {
            // First column is all zero
            x=0;
            ++i;

            for(x=1; x<nx; x++)
            {
                Vec3d pos(x, y, z);
                //pos+=gridmin ;
                if (((P1+i)->data>isoval)^((P1+i-dx)->data>isoval))
                {
                    (P1+i)->p[0] = addPoint(points.wref(), 0, pos,gridmin, (P1+i)->data,(P1+i-dx)->data, mstep, isoval);
                }
                if (((P1+i)->data>isoval)^((P1+i-dy)->data>isoval))
                {
                    (P1+i)->p[1] = addPoint(points.wref(), 1, pos,gridmin,(P1+i)->data,(P1+i-dy)->data, mstep, isoval);
                }
                if (((P1+i)->data>isoval)^((P0+i)->data>isoval))
                {
                    (P1+i)->p[2] = addPoint(points.wref(), 2, pos,gridmin,(P1+i)->data,(P0+i)->data, mstep, isoval);
                }

                // All points should now be created

                if ((P0+i-dx-dy)->data > isoval) mk = 1;
                else mk=0;
                if ((P0+i   -dy)->data > isoval) mk|= 2;
                if ((P0+i      )->data > isoval) mk|= 4;
                if ((P0+i-dx   )->data > isoval) mk|= 8;
                if ((P1+i-dx-dy)->data > isoval) mk|= 16;
                if ((P1+i   -dy)->data > isoval) mk|= 32;
                if ((P1+i      )->data > isoval) mk|= 64;
                if ((P1+i-dx   )->data > isoval) mk|= 128;


                tri=sofa::helper::MarchingCubeTriTable[mk];
                while (*tri>=0)
                {
                    typename std::vector<CubeData>::iterator b = base+i;
                    if (addFace(triangles.wref(),
                                (b+edgecube[tri[0]])->p[edgepts[tri[0]]],
                                (b+edgecube[tri[1]])->p[edgepts[tri[1]]],
                                (b+edgecube[tri[2]])->p[edgepts[tri[2]]], points.size())<0)
                    {
                        /*std::cout << "  mk=0x"<<std::hex<<mk<<std::dec<<" p1="<<tri[0]<<" p2="<<tri[1]<<" p3="<<tri[2]<<sendl;
                        for (int e=0; e<12; e++)
                            std::cout << "  e"<<e<<"="<<(b+edgecube[e])->p[edgepts[e]];
                        std::cout<<sendl;
                        */
                    }
                    tri+=3;
                }
                ++i;
            }
        }
    }
}

void SurfaceMeshGenerationFromScalarField::newPlane()
{
    CubeData c;
    c.p[0] = -1;
    c.p[1] = -1;
    c.p[2] = -1;
    c.data = 0;
    typename std::vector<CubeData>::iterator P = P0;
    P0 = P1;
    P1 = P;
    int n = planes.getValue().size()/2;
    for (int i=0; i<n; ++i,++P)
        *P = c;
}

} /// _surfacemeshgenerationfromscalarfield_
} /// namespace mapping
} /// namespace component
} /// namespace sofa
