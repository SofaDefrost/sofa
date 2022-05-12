#pragma once

#include <SofaImplicitField/config.h>
#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/Ray.h>

#include <sofa/core/DataEngine.h>
#include <SofaImplicitField/components/geometry/DisplacementField.h>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/core/behavior/MechanicalState.h>

#include <SofaOpenglVisual/OglShader.h>
#include <sofa/gl/FrameBufferObject.h>

namespace sofa::component::engine
{

using sofa::core::objectmodel::BaseLink;
using sofa::core::objectmodel::SingleLink;
using sofaimplicitfield::DisplacementField;
using sofa::core::DataEngine;
using sofa::type::Vec3;
using sofa::type::Vec2i;

// Make Template with <FieldIn, FieldOut> with Variations:
//      <DisplacementField, DisplacementField>
//      <DisplacementField, ScalarField>
//      <ScalarField, DisplacementField>
// ?
class SOFA_SOFAIMPLICITFIELD_API ImageBasedVolumeEngineOgl : public DataEngine
{
public:
    SOFA_CLASS(ImageBasedVolumeEngineOgl, DataEngine);
    
    void init() override;
    void reinit() override;
    void doUpdate() override;

    // Inputs:
    SingleLink<ImageBasedVolumeEngineOgl, DisplacementField, BaseLink::FLAG_STRONGLINK> l_field_one;
    SingleLink<ImageBasedVolumeEngineOgl, sofa::component::visualmodel::OglShader, BaseLink::FLAG_STOREPATH> l_shader_one;
    SingleLink<ImageBasedVolumeEngineOgl, DisplacementField, BaseLink::FLAG_STRONGLINK> l_field_two;
    SingleLink<ImageBasedVolumeEngineOgl, sofa::component::visualmodel::OglShader, BaseLink::FLAG_STOREPATH> l_shader_two;
    Data<Vec2i> d_resolution;
    Data<double> d_epsilon;
    // Outputs:
    Data<sofa::helper::vector<Vec3>> d_intersections;
    Data<double> d_volume;
    Data<sofa::helper::vector<Vec3>> d_volume_gradients_one;
    Data<sofa::helper::vector<Vec3>> d_volume_gradients_two;

protected:
    ImageBasedVolumeEngineOgl();
    ~ImageBasedVolumeEngineOgl() override {}

    class InternalData;
    std::unique_ptr<InternalData> data_one;
    std::unique_ptr<InternalData> data_two;
};

}