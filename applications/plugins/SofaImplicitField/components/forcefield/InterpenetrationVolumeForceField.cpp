#include <sofa/core/ObjectFactory.h>

#include "InterpenetrationVolumeForceField.h"
#include "InterpenetrationVolumeForceField.inl"

#include <sofa/defaulttype/VecTypes.h>

namespace sofa::component::forcefield
{

using namespace sofa::defaulttype;

/// Register in the Factory
static int InterpenetrationVolumeForceFieldClass = core::RegisterObject("Repultion ForceField reliant on the interpenetration volume of two DisplacementFields.")
    .add< InterpenetrationVolumeForceField<Vec3Types> >();
template class SOFA_SOFAIMPLICITFIELD_API InterpenetrationVolumeForceField<Vec3Types>;

} // namespace sofa:component:forcefield