#pragma once

#include <SofaImplicitField/config.h>
#include <SofaImplicitField/components/geometry/DisplacementField.h>

#include <sofa/core/objectmodel/DataFileName.h>
#include <sofa/core/objectmodel/Data.h>

#include <sofa/core/behavior/ForceField.h>

#include <sofa/core/behavior/MechanicalState.h>

#include <sofa/defaulttype/VecTypes.h>

namespace sofa
{

namespace component
{

namespace forcefield
{

using sofa::core::objectmodel::BaseLink;
using sofa::core::objectmodel::SingleLink;
using sofaimplicitfield::DisplacementField;
using sofa::type::Vec3;

template<class DataTypes>
class InterpenetrationVolumeForceField : public core::behavior::ForceField<DataTypes>
{
public:

    SOFA_CLASS(SOFA_TEMPLATE(InterpenetrationVolumeForceField, DataTypes), SOFA_TEMPLATE(core::behavior::ForceField, DataTypes));

    typedef core::behavior::ForceField<DataTypes> Inherit;
    typedef typename DataTypes::Real Real;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Deriv Deriv;
    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef Data<VecCoord> DataVecCoord;
    typedef Data<VecDeriv> DataVecDeriv;

    SingleLink<InterpenetrationVolumeForceField, DisplacementField, BaseLink::FLAG_STRONGLINK> l_field;
    Data<double> d_volume;
    Data<sofa::helper::vector<Vec3>> d_volume_gradients;
    Data<double> d_multiplication_factor_k;

    void init() override;
    void addForce(const core::MechanicalParams* params, DataVecDeriv& f, const DataVecCoord& x, const DataVecDeriv& v) override;
    SReal getPotentialEnergy(const core::MechanicalParams* params, const DataVecCoord& x) const override;
    void draw(const core::visual::VisualParams* vparams) override;
    // Implemented but left empty:
    void addDForce(const core::MechanicalParams* mparams, DataVecDeriv& d_df , const DataVecDeriv& d_dx) override;
    /*
    void addKToMatrix(sofa::defaulttype::BaseMatrix *mat, SReal k, unsigned int &offset) override;
    */

protected:
    InterpenetrationVolumeForceField();
    ~InterpenetrationVolumeForceField() override {}
};

} // namespace forcefield

} // namespace component

} // namespace sofa