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
#include <sofa/helper/Factory.inl>
#include <sofa/core/collision/Contact.h>
using sofa::core::collision::Contact ;

#include <SofaBaseCollision/CubeModel.h>
#include <SofaBaseCollision/SphereModel.h>

#include <SofaMeshCollision/RigidContactMapper.inl>
#include <SofaMeshCollision/PointModel.h>
#include <SofaMeshCollision/LineModel.h>
#include <SofaMeshCollision/TriangleModel.h>

#include <SofaMeshCollision/BarycentricContactMapper.inl>
#include <SofaMeshCollision/IdentityContactMapper.h>
#include <SofaBaseCollision/BaseContactMapper.h>
#include <SofaConstraint/FrictionContact.inl>
#include <SofaConstraint/BarycentricDistanceLMConstraintContact.inl>
#include <SofaMiscCollision/BarycentricStickContact.inl>
#include <SofaMeshCollision/BarycentricPenalityContact.inl>
#include <SofaMeshCollision/BarycentricContactMapper.h>
#include <SofaMeshCollision/BarycentricContactMapper.inl>

#include "components/collision/ScalarFieldCollisionModel.h"

namespace sofa::component::collision
{

/// Mapper for FFDDistanceGridCollisionModel
template <class DataTypes>
class ContactMapper<ScalarFieldCollisionModel, DataTypes> : public RigidContactMapper<ScalarFieldCollisionModel, DataTypes>
{
public:
    typedef typename DataTypes::Real Real;
    typedef typename DataTypes::Coord Coord;
};
} /// namespace sofa::component::collision

namespace
{
using sofa::component::collision::ScalarFieldCollisionModel;
using sofa::component::collision::BarycentricPenalityContact;
using sofa::component::collision::BarycentricStickContact;
using sofa::component::collision::BarycentricContactMapper;
using sofa::component::collision::BarycentricDistanceLMConstraintContact;
using sofa::component::collision::FrictionContact;
using sofa::component::collision::ContactMapper;
using sofa::component::collision::PointCollisionModel;
using sofa::component::collision::LineCollisionModel;

//using sofa::component::collision::SphereCollisionModel;
//using sofa::component::collision::TriangleCollisionModel;

///Creator<Contact::Factory, BarycentricPenalityContact<ScalarFieldCollisionModel, PointCollisionModel<sofa::defaulttype::Vec3Types>> > DistanceGridDistanceGridContactClass("default", true);
Creator<Contact::Factory, BarycentricPenalityContact<ScalarFieldCollisionModel, LineCollisionModel<sofa::defaulttype::Vec3Types>> > DistanceGridDistanceGridContactClass("default", true);

///////////////////////////////// BARYCENTRICSTICK /////////////////////////////////////////////////
///Creator<Contact::Factory, BarycentricStickContact<ScalarFieldCollisionModel, PointCollisionModel<sofa::defaulttype::Vec3Types>> > DistanceGridPointStickContactClass("stick", true);
Creator<Contact::Factory, BarycentricStickContact<ScalarFieldCollisionModel, LineCollisionModel<sofa::defaulttype::Vec3Types>> > DistanceGridPointStickContactClass("stick", true);
//Creator<Contact::Factory, BarycentricStickContact<ScalarFieldCollisionModel, SphereCollisionModel<sofa::defaulttype::Vec3Types>> > DistanceGridSphereStickContactClass("stick", true);
//Creator<Contact::Factory, BarycentricStickContact<ScalarFieldCollisionModel, TriangleCollisionModel<sofa::defaulttype::Vec3Types>> > DistanceGridTriangleStickContactClass("stick", true);

/////////////////////////////////////// FRICTION ///////////////////////////////////////////////////
Creator<Contact::Factory, FrictionContact<ScalarFieldCollisionModel, PointCollisionModel<sofa::defaulttype::Vec3Types>> > DistanceGridPointFrictionContactClass("FrictionContact", true);
//Creator<Contact::Factory, FrictionContact<ScalarFieldCollisionModel, SphereCollisionModel<sofa::defaulttype::Vec3Types>> > DistanceGridSphereFrictionContactClass("FrictionContact", true);
//Creator<Contact::Factory, FrictionContact<ScalarFieldCollisionModel, TriangleCollisionModel<sofa::defaulttype::Vec3Types>> > DistanceGridTriangleFrictionContactClass("FrictionContact", true);

/////////////////////////////////////// BarycentricDistanceLMConstraint ///////////////////////////////////
///Creator<Contact::Factory, BarycentricDistanceLMConstraintContact<ScalarFieldCollisionModel, PointCollisionModel<sofa::defaulttype::Vec3Types>> > DistanceGridPointDistanceLMConstraintContactClass("distanceLMConstraint",true);
Creator<Contact::Factory, BarycentricDistanceLMConstraintContact<ScalarFieldCollisionModel, LineCollisionModel<sofa::defaulttype::Vec3Types>> > DistanceGridPointDistanceLMConstraintContactClass("distanceLMConstraint",true);
//Creator<Contact::Factory, BarycentricDistanceLMConstraintContact<ScalarFieldCollisionModel, SphereCollisionModel<sofa::defaulttype::Vec3Types>> > DistanceGridSphereDistanceLMConstraintContactClass("distanceLMConstraint",true);
//Creator<Contact::Factory, BarycentricDistanceLMConstraintContact<ScalarFieldCollisionModel, TriangleCollisionModel<sofa::defaulttype::Vec3Types>> > DistanceGridTriangleDistanceLMConstraintContactClass("distanceLMConstraint",true);

}
