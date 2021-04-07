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
#pragma once

#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/defaulttype/Vec3Types.h>
#include <sofa/core/CollisionModel.h>
#include <SofaImplicitField/components/geometry/ScalarField.h>
#include <sofa/core/behavior/MechanicalState.h>

namespace sofa::component::collision
{

/// This is the begining of a private namespace (indicated by the _ symbole around the name) .
namespace _scalarfieldcollisionmodel_
{
using sofa::component::geometry::ScalarField;
using sofa::core::TCollisionElementIterator;
using sofa::core::CollisionModel;
using sofa::core::objectmodel::SingleLink;
using sofa::core::behavior::MechanicalState;
using sofa::defaulttype::Rigid3Types;
using sofa::defaulttype::Vec3Types;


class ScalarFieldCollisionModel;
class SOFA_SOFAIMPLICITFIELD_API ScalarFieldCollisionElement : public TCollisionElementIterator<ScalarFieldCollisionModel>
{
public:
    explicit ScalarFieldCollisionElement(const core::CollisionElementIterator& i);

    ScalarFieldCollisionElement(ScalarFieldCollisionModel* model, int index);
    ScalarField* getScalarField();
};

class SOFA_SOFAIMPLICITFIELD_API ScalarFieldCollisionModel : public CollisionModel
{
public:
    SOFA_CLASS(ScalarFieldCollisionModel, CollisionModel);
    typedef ScalarFieldCollisionElement Element;
    typedef Rigid3Types InDataTypes;
    typedef Vec3Types DataTypes;

    SingleLink<ScalarFieldCollisionModel, ScalarField,
               BaseLink::FLAG_STOREPATH | BaseLink::FLAG_STRONGLINK> l_field;

    MechanicalState<InDataTypes>* getMechanicalState() { return rigid; }
    MechanicalState<Rigid3Types>* rigid;

    void init() override;
    void draw(const sofa::core::visual::VisualParams*) override;
    ScalarField* getScalarField() const { return l_field.get(); }

protected:
    ScalarFieldCollisionModel();

    void computeBoundingTree(int maxDepth=0) override ;
};


} /// namespace _scalarfieldcollisionmodel_

using _scalarfieldcollisionmodel_::ScalarFieldCollisionModel;
using _scalarfieldcollisionmodel_::ScalarFieldCollisionElement;

} /// namespace sofa::component::collision

