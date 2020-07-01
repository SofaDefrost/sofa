
#include <SofaImplicitField/components/collision/ScalarFieldIntersectionLine.h>
#include <SofaBaseCollision/DiscreteIntersection.h>
#include <sofa/core/collision/IntersectorFactory.h>
#include <sofa/core/collision/Intersection.inl>
#include <sofa/helper/logging/Messaging.h>
#include <iostream>
using namespace std;
namespace sofa::component::collision::_scalarfieldintersectionline_
{
using sofa::core::collision::IntersectorCreator;
IntersectorCreator<DiscreteIntersection, ScalarFieldIntersectionLine> ScalarGridIntersector("ScalarFieldIntersection");

ScalarFieldIntersectionLine::ScalarFieldIntersectionLine(DiscreteIntersection* object) :
    intersection(object)
{
    intersection->intersectors.add<ScalarFieldCollisionModel, LineCollisionModel<sofa::defaulttype::Vec3Types>,
                                   ScalarFieldIntersectionLine>(this);
}

bool ScalarFieldIntersectionLine::testIntersection(ScalarFieldCollisionElement&, Line&)
{
    return true;
}



//double ScalarFieldIntersectionLine::Interpolation(double a , double b , double lambda){

//    return a + lambda*(b-a);
//}
Vec3d ScalarFieldIntersectionLine::Interpolation_point(Vec3d& p1, Vec3d& p2, double lambda){

    return p1 +(p2-p1)*lambda ;
}

int ScalarFieldIntersectionLine::computeIntersection(ScalarFieldCollisionElement& collisionElement,
                                                 Line& line, OutputVector* contacts)
{

//    ///cout << "line";
//    msg_warning("cc") << "line ";
//    /// Get the scalar field out of the collision element
       sofa::component::geometry::ScalarField* field = collisionElement.getScalarField();

//    /// In case there is no scalar field associated to this object we can just leave.
//    if(!field)
//        return 0;

//    /// retrieving the 2 points of the Line line.
//    sofa::defaulttype::Vec3 p1 = line.p1();
//    sofa::defaulttype::Vec3 p2 = line.p2();
//    ///sofa::defaulttype::Vec3 p_previous = line.p(i);
//    ///sofa::defaulttype::Vec3 p_next = line.p(i+1);


//    /// Evaluate the field at points p1 & p2.
//    SReal d1 = field->getValue(p1);
//    SReal d2 = field->getValue(p2);

//    /// Compute the common proximity distance.
//    const double d0 = line.getProximity() + intersection->getContactDistance();
//    const SReal margin = SReal{0.001} + SReal{d0};

//    /// If the implicit field's value is bigger than zero + the margin then we exit without collision
//    /// Case wwhere we have no collision for the 2 points
//    if (d1 >= margin && d2 >= margin)
//        return 0;

//    /// Here we are in collision so we will compute the gradient to estimate the
//    /// direction to push out the penetrating point.
//    /// Case where we have p1 colliding
//    else if (d1 < margin ){
//        sofa::defaulttype::Vec3 grad1 = field->getGradient(p1);

//        /// Generates the contact response and add it to the output set of contacts.
//        contacts->resize(contacts->size()+1);
//        core::collision::DetectionOutput *detection = &*(contacts->end()-1);
//        detection->normal = grad1;       ///< normal in global space from p1's surface
//        detection->value = (d1- d0);     ///< penetrating distance intensity
//        detection->elem.first = collisionElement;     ///< first colliding element
//        detection->elem.second = line;               ///< second colliding element
//        detection->id = line.getIndex();
//        detection->point[0] = p1 - grad1 * detection->value;
//        detection->point[1] = p1;

//        return 1; ///< returns the number of detected collision point.
//    }
//    ///Case where we have p2 colliding
//    else if(d2 < margin){
//    sofa::defaulttype::Vec3 grad2 = field->getGradient(p2);

//    /// Generates the contact response and add it to the output set of contacts.
//    contacts->resize(contacts->size()+1);
//    core::collision::DetectionOutput *detection = &*(contacts->end()-1);
//    detection->normal = grad2;       ///< normal in global space from p1's surface
//    detection->value = (d2- d0);     ///< penetrating distance intensity
//    detection->elem.first = collisionElement;     ///< first colliding element
//    detection->elem.second = line;               ///< second colliding element
//    detection->id = line.getIndex();
//    detection->point[0] = p2 - grad2 * detection->value;
//    detection->point[1] = p2;
//    return 1; ///< returns the number of detected collision point.
//}






        for(int i=2;i<10;i++){
          double alpha = i/9.0 ;
          sofa::defaulttype::Vec3 p0 = line.p1();
          sofa::defaulttype::Vec3 p1 = line.p2();
          Vec3d p_i = Interpolation_point(p0,p1,alpha);
          SReal d_i = field->getValue(p_i);
          const double di = line.getProximity() + intersection->getContactDistance();
          const SReal margin_i = SReal{0.001} + SReal{di};
          if (di >= margin_i)
              return 0;
          else if (di < margin_i ){
              sofa::defaulttype::Vec3 grad_i = field->getGradient(p_i);

              /// Generates the contact response and add it to the output set of contacts.
              contacts->resize(contacts->size()+1);
              core::collision::DetectionOutput *detection = &*(contacts->end()-1);
              detection->normal = grad_i;       ///< normal in global space from p1's surface
              detection->value = (d_i- di);     ///< penetrating distance intensity
              detection->elem.first = collisionElement;     ///< first colliding element
              detection->elem.second = line;               ///< second colliding element
              detection->id = line.getIndex();
              detection->point[0] = p_i - grad_i * detection->value;
              detection->point[1] = p_i;

              return i; ///< returns the number of detected collision point. check this
          }


        }

        return 0 ;
    }

};

