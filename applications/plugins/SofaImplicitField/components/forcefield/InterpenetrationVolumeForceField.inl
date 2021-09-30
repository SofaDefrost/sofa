#pragma once

#include <sofa/core/visual/VisualParams.h>
#include <sofa/core/behavior/MultiMatrixAccessor.h>
#include "InterpenetrationVolumeForceField.h"


namespace sofa::component::forcefield
{

using sofa::helper::getReadAccessor;
using sofa::helper::getWriteAccessor;

template<class DataTypes>
InterpenetrationVolumeForceField<DataTypes>::InterpenetrationVolumeForceField():
    l_field(initLink("field", "The targeted DisplacementField.")),
    d_volume(initData(&d_volume, "volume", "The evaluated interpenetration volume.")),
    d_volume_gradients(initData(&d_volume_gradients, "volume_gradients", "The evaluated interpenetration volume' gradients w.r.t. the DOFs of the DisplacementField.")),
    d_multiplication_factor_k(initData(&d_multiplication_factor_k, 1.0, "k", "A positive multiplication factor."))
{
    //TODO ?
    msg_warning() << "InterpenetrationVolumeForceField<DataTypes>::InterpenetrationVolumeForceField()";
}

template<class DataTypes>
void InterpenetrationVolumeForceField<DataTypes>::init()
{
    msg_warning() << "InterpenetrationVolumeForceField<DataTypes>::init()";

    this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);

    if (l_field.empty())
    {
        msg_warning() << "ERROR";
        return;
    }
    // Get the MechanicalState associated with the passed field.
    core::behavior::MechanicalState<DataTypes>* _mechaState = l_field->l_dofs.get();
    //this->getContext()->addObject(_mechaState);

    // Call the init() from ForceField.
    this->mstate.set(_mechaState); //this->mstate.set(dynamic_cast< core::behavior::MechanicalState<DataTypes>* >(l_field->l_dofs.get()));
    Inherit::init();

    // add to tracker
    this->trackInternalData(d_volume);
    this->trackInternalData(d_volume_gradients);
    this->trackInternalData(d_multiplication_factor_k);

    // if all init passes, component is valid
    this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Valid);
}

template<class DataTypes>
void InterpenetrationVolumeForceField<DataTypes>::addForce(const core::MechanicalParams* params, DataVecDeriv& f, const DataVecCoord& x, const DataVecDeriv& v)
{
    msg_warning() << "InterpenetrationVolumeForceField<DataTypes>::addForce(...)";

    SOFA_UNUSED(params);
    SOFA_UNUSED(x);
    SOFA_UNUSED(v);

    auto volume = getReadAccessor(d_volume);
    auto volume_gradients = getReadAccessor(d_volume_gradients);
    auto k = getReadAccessor(d_multiplication_factor_k);

    auto _f = getWriteAccessor(f);

    for (unsigned int i=0; i<volume_gradients.size(); i++)
    {
        _f[i] += k * volume * volume_gradients[i];
    }
}

template<class DataTypes>
SReal InterpenetrationVolumeForceField<DataTypes>::getPotentialEnergy(const core::MechanicalParams* params, const DataVecCoord& x) const
{
    SOFA_UNUSED(params);
    SOFA_UNUSED(x);

    double volume = getReadAccessor(d_volume);
    double k = getReadAccessor(d_multiplication_factor_k);

    return 0.5 * k * pow(volume,2);
}

template<class DataTypes>
void InterpenetrationVolumeForceField<DataTypes>::addDForce(const core::MechanicalParams* mparams, DataVecDeriv& d_df , const DataVecDeriv& d_dx)
{
    // Derivative of a constant force is null, no need to compute addKToMatrix nor addDForce
    SOFA_UNUSED(mparams);
    SOFA_UNUSED(d_df);
    SOFA_UNUSED(d_dx);
    // mparams->setKFactorUsed(true);
}
/*
template<class DataTypes>
void InterpenetrationVolumeForceField<DataTypes>::addKToMatrix(sofa::defaulttype::BaseMatrix * mat, SReal k, unsigned int & offset)
{
    // Derivative of a constant force is null, no need to compute addKToMatrix nor addDForce
    SOFA_UNUSED(mat);
    SOFA_UNUSED(k);
    SOFA_UNUSED(offset);
}
*/
template<class DataTypes>
void InterpenetrationVolumeForceField<DataTypes>::draw(const core::visual::VisualParams* vparams)
{
    if (!vparams->displayFlags().getShowForceFields()) return;
    
    // TODO
}


} // namespace sofa:component:forcefield