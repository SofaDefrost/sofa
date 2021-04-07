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

#include <SofaImplicitField/components/collision/ScalarFieldIntersection.h>
#include <SofaBaseCollision/DiscreteIntersection.h>
#include <sofa/core/collision/IntersectorFactory.h>
#include <sofa/core/collision/Intersection.inl>
namespace sofa::component::collision::_scalarfieldintersection_
{
using sofa::core::collision::IntersectorCreator;
IntersectorCreator<DiscreteIntersection, ScalarFieldIntersection> ScalarGridIntersector("ScalarFieldIntersection");

ScalarFieldIntersection::ScalarFieldIntersection(DiscreteIntersection* object) :
    intersection(object)
{
    intersection->intersectors.add<ScalarFieldCollisionModel, PointCollisionModel<sofa::defaulttype::Vec3Types>,
                                   ScalarFieldIntersection>(this);
}

bool ScalarFieldIntersection::testIntersection(ScalarFieldCollisionElement&, Point&)
{
    return true;
}

int ScalarFieldIntersection::computeIntersection(ScalarFieldCollisionElement& collisionElement,
                                                 Point& point, OutputVector* contacts)
{
    /// Get the scalar field out of the collision element
    sofa::component::geometry::ScalarField* field = collisionElement.getScalarField();

    /// In case there is no scalar field associated to this object we can just leave.
    if(!field)
        return 0;

    /// The the evaluation point.
    sofa::defaulttype::Vec3 p1 = point.p();

    /// Evaluate the field at this p1 point.
    SReal d = field->getValue(p1);

    /// Compute the common proximity distance.
    const double d0 = point.getProximity() + intersection->getContactDistance();
    const SReal margin = SReal{0.001} + SReal{d0};

    /// If the implicit field's value is bigger than zero + the margin then we exit without collision
    if (d >= margin)
        return 0;

    /// Here we are in collision so we will compute the gradient to estimate the
    /// direction to push out the penetrating point.
    sofa::defaulttype::Vec3 grad = field->getGradient(p1);

    /// Generates the contact response and add it to the output set of contacts.
    contacts->resize(contacts->size()+1);
    core::collision::DetectionOutput *detection = &*(contacts->end()-1);
    detection->normal = grad;       ///< normal in global space from p1's surface
    detection->value = (d- d0);     ///< penetrating distance intensity
    detection->elem.first = collisionElement;     ///< first colliding element
    detection->elem.second = point;               ///< second colliding element
    detection->id = point.getIndex();
    detection->point[0] = p1 - grad * detection->value;
    detection->point[1] = p1;

    return 1; ///< returns the number of detected collision point.
}

};
