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
#ifndef SOFA_COMPONENT_COLLISION_MORFRICTIONCONTACT_H
#define SOFA_COMPONENT_COLLISION_MORFRICTIONCONTACT_H
#include "config.h"

#include <SofaConstraint/FrictionContact.h>
#include <sofa/core/collision/Contact.h>
#include <sofa/core/collision/Intersection.h>
#include <SofaBaseMechanics/BarycentricMapping.h>
#include <SofaConstraint/UnilateralInteractionConstraint.h>
#include <SofaConstraint/MORUnilateralInteractionConstraint.h>
#include <sofa/helper/Factory.h>
#include <SofaBaseCollision/BaseContactMapper.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/BaseMapping.h>

#include <SofaConstraint/ContactIdentifier.h>

namespace sofa
{

namespace component
{

namespace collision
{



template <class TCollisionModel1, class TCollisionModel2, class ResponseDataTypes = sofa::defaulttype::Vec3Types >
class MORFrictionContact : public FrictionContact<TCollisionModel1, TCollisionModel2, ResponseDataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE3(MORFrictionContact, TCollisionModel1, TCollisionModel2, ResponseDataTypes), SOFA_TEMPLATE3(FrictionContact,TCollisionModel1, TCollisionModel2, ResponseDataTypes));
    typedef TCollisionModel1 CollisionModel1;
    typedef TCollisionModel2 CollisionModel2;
    typedef core::collision::Intersection Intersection;
    typedef typename TCollisionModel1::DataTypes::CPos TVec1;
    typedef typename TCollisionModel1::DataTypes::CPos TVec2;
    typedef sofa::defaulttype::StdVectorTypes<TVec1, TVec2, typename TCollisionModel1::DataTypes::Real > DataTypes1;
    typedef sofa::defaulttype::StdVectorTypes<TVec1, TVec2, typename TCollisionModel1::DataTypes::Real > DataTypes2;

    typedef core::behavior::MechanicalState<DataTypes1> MechanicalState1;
    typedef core::behavior::MechanicalState<DataTypes2> MechanicalState2;
    typedef typename CollisionModel1::Element CollisionElement1;
    typedef typename CollisionModel2::Element CollisionElement2;
    typedef core::collision::DetectionOutputVector OutputVector;
    typedef core::collision::TDetectionOutputVector<CollisionModel1,CollisionModel2> TOutputVector;

protected:
    using FrictionContact<TCollisionModel1, TCollisionModel2, ResponseDataTypes>::model1;
    using FrictionContact<TCollisionModel1, TCollisionModel2, ResponseDataTypes>::model2;
    using FrictionContact<TCollisionModel1, TCollisionModel2, ResponseDataTypes>::intersectionMethod;
    using FrictionContact<TCollisionModel1, TCollisionModel2, ResponseDataTypes>::selfCollision; ///< true if model1==model2 (in this case, only mapper1 is used)
    using FrictionContact<TCollisionModel1, TCollisionModel2, ResponseDataTypes>::mapper1;
    using FrictionContact<TCollisionModel1, TCollisionModel2, ResponseDataTypes>::mapper2;

    using FrictionContact<TCollisionModel1, TCollisionModel2, ResponseDataTypes>::m_constraint;

    using FrictionContact<TCollisionModel1, TCollisionModel2, ResponseDataTypes>::parent;

    using FrictionContact<TCollisionModel1, TCollisionModel2, ResponseDataTypes>::mu; ///< friction coefficient (0 for frictionless contacts)
    using FrictionContact<TCollisionModel1, TCollisionModel2, ResponseDataTypes>::tol; ///< tolerance for the constraints resolution (0 for default tolerance)
    using FrictionContact<TCollisionModel1, TCollisionModel2, ResponseDataTypes>::contacts;
    using FrictionContact<TCollisionModel1, TCollisionModel2, ResponseDataTypes>::mappedContacts;

    void activateMappers() override;
    MORFrictionContact();
    MORFrictionContact(CollisionModel1* model1_, CollisionModel2* model2_, Intersection* intersectionMethod_);
    ~MORFrictionContact() override;

//    void setInteractionTags(MechanicalState1* mstate1, MechanicalState2* mstate2);

//    FrictionContact();
//    FrictionContact(CollisionModel1* model1, CollisionModel2* model2, Intersection* intersectionMethod);

//    ~FrictionContact() override;
public:
    void cleanup() override;

//    std::pair<core::CollisionModel*,core::CollisionModel*> getCollisionModels() override { return std::make_pair(model1,model2); }

//    void setDetectionOutputs(OutputVector* outputs) override;

//    void createResponse(core::objectmodel::BaseContext* group) override;

//    void removeResponse() override;
};

} // collision

} // component

} // sofa

#endif // SOFA_COMPONENT_COLLISION_MORFRICTIONCONTACT_H
