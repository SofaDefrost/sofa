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
#define SOFA_COMPONENT_COLLISION_DISTANCEGRIDCOLLISIONMODEL_CPP

#include<SofaImplicitField/components/collision/ScalarFieldCollisionModel.h>
#include <SofaBaseCollision/CubeModel.h>
#include<sofa/core/ObjectFactory.h>
using sofa::core::RegisterObject;

namespace sofa::component::collision::_scalarfieldcollisionmodel_
{

ScalarFieldCollisionElement::ScalarFieldCollisionElement(ScalarFieldCollisionModel* model, int index) :
    TCollisionElementIterator<ScalarFieldCollisionModel> (model, index)
{}

ScalarFieldCollisionElement::ScalarFieldCollisionElement(const core::CollisionElementIterator& i)
    : TCollisionElementIterator<ScalarFieldCollisionModel>(static_cast<ScalarFieldCollisionModel*>(i.getCollisionModel()), i.getIndex())
{
}

ScalarField* ScalarFieldCollisionElement::getScalarField()
{
    return dynamic_cast<ScalarField*>(model->getScalarField());
}

ScalarFieldCollisionModel::ScalarFieldCollisionModel() :
    l_field(initLink("field","An object implementing the ScalarField interface."))
{

}

void ScalarFieldCollisionModel::init()
{
    CollisionModel::init();
    rigid = dynamic_cast< MechanicalState<Rigid3Types>* > (getContext()->getMechanicalState());
    if(rigid == nullptr)
        msg_error() << "Missing Mechanical State" ;
}

void ScalarFieldCollisionModel::computeBoundingTree(int maxDepth)
{
    CubeCollisionModel* cubeModel = this->createPrevious<CubeCollisionModel>();

    if (!cubeModel->empty()) return; // No need to recompute BBox if immobile

    cubeModel->resize(1);
    cubeModel->setParentOf(0, {-6.0,-6.0,-6.0}, {6.0,6.0,6.0}); // define the bounding box of the current element
    cubeModel->computeBoundingTree(maxDepth);
}

void ScalarFieldCollisionModel::draw(const sofa::core::visual::VisualParams* vparms)
{
    if (getPrevious())
        getPrevious()->draw(vparms);
}

/// Register the object in the Factory
static int ScalarFieldCollisionModelClass = RegisterObject("A collision model for implicit fields.")
        .add< ScalarFieldCollisionModel >()
        ;

} /// sofa::component::collision
