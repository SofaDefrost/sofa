#pragma once

#include <SofaImplicitField/config.h>
#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/defaulttype/VecTypes.h>
#include <SofaImplicitField/components/geometry/DisplacementField.h>
#include <SofaImplicitField/components/engine/RayMarchingEngine.h>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/core/behavior/MechanicalState.h>

namespace sofa::component::engine
{

using sofa::core::objectmodel::BaseLink;
using sofa::core::objectmodel::SingleLink;
using sofaimplicitfield::DisplacementField;
using sofa::defaulttype::Vec3;
using sofa::defaulttype::Vec2i;

class SOFA_SOFAIMPLICITFIELD_API ImageBasedVolumeEngine : public RayMarchingEngine
{
public:
    SOFA_CLASS(ImageBasedVolumeEngine, RayMarchingEngine);

    void doUpdate() override;

    // Inputs:
    SingleLink<ImageBasedVolumeEngine, DisplacementField, BaseLink::FLAG_STRONGLINK> l_field_one;
    SingleLink<ImageBasedVolumeEngine, DisplacementField, BaseLink::FLAG_STRONGLINK> l_field_two;
    Data<Vec2i> d_resolution;
    Data<double> d_epsilon;
    // Outputs:
    Data<sofa::helper::vector<Vec3>> d_intersections;
    Data<double> d_volume;
    Data<sofa::helper::vector<Vec3>> d_volume_gradients_one;
    Data<sofa::helper::vector<Vec3>> d_volume_gradients_two;

protected:
    ImageBasedVolumeEngine();
    ~ImageBasedVolumeEngine() override {}

    struct Hit
    {
        bool found = false;
        Vec3 pos;
        Vec3 normal;
        double distance;
        bool surface_id;
        int domain;
    };

    Hit sphereTracing(const sofa::defaulttype::Ray& r, const double eps, const double max_depth);
};

}