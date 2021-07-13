#include <sofa/core/ObjectFactory.h>
#include <sofa/core/objectmodel/Data.h>
#include "ImageBasedVolumeEngine.h"

namespace sofa::component::engine
{
using sofa::defaulttype::Ray;
using sofa::helper::getReadAccessor;
using sofa::helper::getWriteAccessor;
using sofaimplicitfield::DisplacementField;

/// Register in the Factory
static int ImageBasedVolumeEngineClass = core::RegisterObject("Mono-volume setting.").add< ImageBasedVolumeEngine >();

ImageBasedVolumeEngine::ImageBasedVolumeEngine():
    // Inputs:   <-    <-    <-    <-    <-     <-     <-   <-    <-    <-   Why are all of these undefined? 
    l_field_one(initLink("field one", "The first colliding scalar field.")),
    l_topology_one(initLink("topology one", "The mesh topology supporting the first scalar field")),
    l_dofs_one(initLink("dofs one", "The nodal values of the second scallar field.")),
    l_field_two(initLink("field two", "The second colliding scalar field.")),
    l_topology_two(initLink("topology two", "The mesh topology supporting the second scalar field")),
    l_dofs_two(initLink("dofs two ", "The nodal values of the second scallar field.")),
    d_resolution(initData(&d_resolution, Vec2i{20,20}, "resolution", "The amount of samples per visual axis.")),
    d_epsilon(initData(&d_epsilon, 0.01, "epsilon", "The tolerance allowed when evaluating an implicit surface.")),
    // Outputs:
    d_intersections(initData(&d_intersections, "intersections", "The intersection points' locations.")),
    d_volume(initData(&d_volume, "volume", "The evaluated interpenetration volume.")),
    d_volume_gradients_one(initData(&d_volume_gradients_one, "volume gradients one", "The evaluated interpenetration volume' gradients w.r.t. the DOFs of the first scallar field.")),
    d_volume_gradients_two(initData(&d_volume_gradients_two, "volume gradients two", "The evaluated interpenetration volume' gradients w.r.t. the DOFs of the second scallar field."))
{
    addOutput(&d_intersections);
    addOutput(&d_volume);
    addOutput(&d_volume_gradients_one);
    addOutput(&d_volume_gradients_two);
}

bool ImageBasedVolumeEngine::sphereTracing(const Ray& r, Vec3& out_vec, bool& out_ind, Vec2i& out_tetra, DisplacementField* field_one, DisplacementField* field_two, const double eps, const double max_depth)
{
    double travelled = 0.0;
    Vec3 current_pos = r.origin();
    const Vec3& dir = r.direction();
    while (travelled <= max_depth)
    {
        // Naïve search for the parent volumetric primitive.
        int t_one = field_one->getDomain(current_pos);
        int t_two = field_two->getDomain(current_pos);
        // evaluate both implicit functions.
        double dist_one = field_one->getValue(current_pos, t_one);
        double dist_two = field_two->getValue(current_pos, t_two);
        double dist = fabs(fmin(dist_one, dist_two));
        if(dist<eps)
        {
            out_vec = current_pos;
            out_ind = (dist_one<=dist_two) ? true : false;
            out_tetra[0] = t_one;
            out_tetra[1] = t_two;
            return true;
        }
        travelled += dist;
        current_pos = current_pos + dir * dist;
    }
    return false;
}

void RayMarchingEngine::doUpdate()
{
    if(l_field_one.empty() || l_field_two.empty())
    {
        return;
    }
    ///// THEN tell everthing is (will be) up to date now
    /// @warning This must be done AFTER updating all inputs
    /// can be done before or after setting up the outputs
    cleanDirty();

    // Initialize accessors.
    // Inputs.
    DisplacementField* field_one = l_field_one.get();
    DisplacementField* field_two = l_field_two.get();
    auto res = getReadAccessor(d_resolution);
    double eps = getReadAccessor(d_epsilon);
    auto bbox_one = getReadAccessor(*l_dofs_one->read(sofa::defaulttype::BoundingBox::f_bbox()));
    auto dof_one = getReadAccessor(*l_dofs_one->read(sofa::core::VecCoordId::position()));
    auto tetra_one = l_topology_one->getTetrahedronArray();
    auto bbox_two = getReadAccessor(*l_dofs_two->read(sofa::defaulttype::BoundingBox::f_bbox()));
    auto dof_two = getReadAccessor(*l_dofs_two->read(sofa::core::VecCoordId::position()));
    auto tetra_two = l_topology_two->getTetrahedronArray();
    // Outputs.
    auto intersections = getWriteAccessor(d_intersections);
    double volume = getWriteAccessor(d_volume);
    auto volume_gradients_one = getWriteAccessor(d_volume_gradients_one);
    auto volume_gradients_two = getWriteAccessor(d_volume_gradients_two);

    // Clean outputs.
    intersections.clear();
    volume = 0.0;
    volume_gradients_one.clear();
    volume_gradients_one.reserve(dof_one.size());
    volume_gradients_two.clear();
    volume_gradients_two.reserve(dof_two.size());

    // Define sphereTracing output containers.
    Vec3 out_vec;
    bool out_ind;
    Vec2i out_tetra;
    sofa::defaulttype::Vec4d barycentric_coordinates;
    sofa::core::topology::BaseMeshTopology::Tetrahedron tetra;

    /* Broad phase */
    if (bbox_one->intersect(bbox_two))
    {
        // Construct the AABB of the (potential) interpenetration volume.
        auto bbox = bbox_one->getIntersection(bbox_two);
        Vec3 bbox_bottom = bbox->minBBox();
        Vec3 bbox_size = fabs(box->maxBBox()-bbox_bottom);
    }
    else
    {
        msg_warning() << "The broad phase determined that there is no interpenetration.";
        return;
    }
    /* Narrow phase */
    // Iterate over the faces of the AABB, pairing the oposite faces together.
    std::vector<sofa::defaulttype::Vec2i> planes {{1,2}, {0,1}, {0,2}};
    for (unsigned int plane_it=0; plane_it<3; plane_it++)
    {
        // Compute area of each pixel.
        double width = bbox_size[planes[plane_it][0]] / res->x();
        double height = bbox_size[planes[plane_it][1]] / res->y();
        double pixel_area = width * height;
        // Define plane iterator.
        Vec3 next_line {0.0, 0.0, 0.0};
        next_line[planes[plane_it][0]] = width;
        Vec3 next_column {0.0, 0.0, 0.0};
        next_column[planes[plane_it][1]] = height;
        // Define viewing direction.
        Vec3 viewing_direction {1.0, 1.0, 1.0};
        viewing_direction[planes[plane_it][0]] = 0.0;
        viewing_direction[planes[plane_it][1]] = 0.0;
        // Define max depth.
        double max_depth = (bbox_size * viewing_direction).sum();
        // Begin ray casting.
        Vec3 current_line = bbox_bottom;
        for (unsigned int i=0; i<res->x(); i++)
        {
            Vec3 current_column = current_line;
            for (unsigned int j=0; j<res->y(); j++)
            {
                // Launch ray.
                Ray ray {current_column, viewing_direction};
                if (sphereTracing(ray, out_vec, out_ind, out_tetra, field_one, field_two, eps, max_depth))
                {
                    // Store interesection point (for display purposes).
                    intersections.push_back(out_vec);
                    // Accumulate Volume.
                    volume -= pixel_area*((out_vec-current_column)*viewing_direction).sum();
                    // Accumulate volume gradients.
                    if (out_ind)
                    {
                        barycentric_coordinates = field_one->getBarycentricCoordinates(out_vec, out_tetra[0]);
                        tetra = tetra_one[out_tetra[0]];
                        for (unsigned int k=0; k<4; k++)
                        {
                            volume_gradients_one[out_tetra[0]] += pixel_area * barycentric_coordinates[k] * dof_one[tetra[k]] * viewing_direction * -1;
                        }
                    }
                    else
                    {
                        barycentric_coordinates = field_two->getBarycentricCoordinates(out_vec, out_tetra[1]);
                        tetra = tetra_two[out_tetra[1]];
                        for (unsigned int k=0; k<4; k++)
                        {
                            volume_gradients_two[out_tetra[1]] += pixel_area * barycentric_coordinates[k] * dof_two[tetra[k]] * viewing_direction * -1;
                        }
                    }                    
                    // Launch ray from the oposite side of the bounding box.
                    ray.updateRay(current_column+max_depth*viewing_direction, -1*viewing_direction);
                    if (sphereTracing(ray, out_vec, out_ind, out_tetra, field_one, field_two, eps, max_depth))
                    {
                        // Store interesection point (for display purposes).
                        intersections.push_back(out_vec);
                        // Accumulate Volume.
                        volume += pixel_area*((out_vec-current_column)*viewing_direction).sum();
                        // Accumulate volume gradients.
                        if (out_ind)
                        {
                            barycentric_coordinates = field_one->getBarycentricCoordinates(out_vec, out_tetra[0]);
                            tetra = tetra_one[out_tetra[0]];
                            for (unsigned int k=0; k<4; k++)
                            {
                                volume_gradients_one[out_tetra[0]] += pixel_area * barycentric_coordinates[k] * dof_one[tetra[k]] * viewing_direction;
                            }
                        }
                        else
                        {
                            barycentric_coordinates = field_two->getBarycentricCoordinates(out_vec, out_tetra[1]);
                            tetra = tetra_two[out_tetra[1]];
                            for (unsigned int k=0; k<4; k++)
                            {
                                volume_gradients_two[out_tetra[1]] += pixel_area * barycentric_coordinates[k] * dof_two[tetra[k]] * viewing_direction;
                            }
                        }
                    }
                }
                current_column += next_column;
            }
            current_line += next_line;
        }
    }
    volume /= 3;
    return;
}

} /// namespace sofa::component::engine