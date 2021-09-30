#include <sofa/core/ObjectFactory.h>
#include "ImageBasedVolumeEngine.h"

#include <time.h> // Print execusion time of each time step!

namespace sofa::component::engine
{
using sofa::defaulttype::Ray;
using sofa::helper::getReadAccessor;
using sofa::helper::getWriteAccessor;
using sofaimplicitfield::DisplacementField;

/// Register in the Factory
static int ImageBasedVolumeEngineClass = core::RegisterObject("Mono-volume setting. Sequential implementation.").add< ImageBasedVolumeEngine >();

ImageBasedVolumeEngine::ImageBasedVolumeEngine():
    // Inputs:
    l_field_one(initLink("field_one", "The first colliding scalar field.")),
    l_field_two(initLink("field_two", "The second colliding scalar field.")),
    d_resolution(initData(&d_resolution, Vec2i{20,20}, "resolution", "The amount of samples per visual axis.")),
    d_epsilon(initData(&d_epsilon, 0.01, "epsilon", "The tolerance allowed when evaluating an implicit surface.")),
    // Outputs:
    d_intersections(initData(&d_intersections, "intersections", "The intersection points' locations.")),
    d_volume(initData(&d_volume, "volume", "The evaluated interpenetration volume.")),
    d_volume_gradients_one(initData(&d_volume_gradients_one, "volume_gradients_one", "The evaluated interpenetration volume' gradients w.r.t. the DOFs of the first scallar field.")),
    d_volume_gradients_two(initData(&d_volume_gradients_two, "volume_gradients_two", "The evaluated interpenetration volume' gradients w.r.t. the DOFs of the second scallar field."))
{
    addOutput(&d_intersections);
    addOutput(&d_volume);
    addOutput(&d_volume_gradients_one);
    addOutput(&d_volume_gradients_two);
}

void ImageBasedVolumeEngine::init()
{
    setDirtyValue();

    
    // Initialize outputs to 0.
    auto dof_one = getReadAccessor(*l_field_one->l_dofs->read(sofa::core::VecCoordId::position()));
    auto dof_two = getReadAccessor(*l_field_two->l_dofs->read(sofa::core::VecCoordId::position()));

    auto volume_gradients_one = getWriteAccessor(d_volume_gradients_one);
    auto volume_gradients_two = getWriteAccessor(d_volume_gradients_two);
    
    volume_gradients_one.clear();
    volume_gradients_one.resize(dof_one.size());
    std::fill(volume_gradients_one.begin(), volume_gradients_one.end(), Vec3{0.0, 0.0, 0.0});
    volume_gradients_two.clear();
    volume_gradients_two.resize(dof_two.size());
    std::fill(volume_gradients_two.begin(), volume_gradients_two.end(), Vec3{0.0, 0.0, 0.0});

    d_volume.beginEdit();
    d_volume.setValue(0.0);
    d_volume.endEdit();

    // add to tracker ?
    this->trackInternalData(d_intersections);
    this->trackInternalData(d_volume);
    this->trackInternalData(d_volume_gradients_one);
    this->trackInternalData(d_volume_gradients_two);

    // cleanDirty();
}

void ImageBasedVolumeEngine::reinit()
{
    setDirtyValue();
    update();
}

ImageBasedVolumeEngine::Hit ImageBasedVolumeEngine::sphereTracing(const sofa::defaulttype::Ray& r, const double eps, const double max_depth)
{
    double travelled = 0.0;
    Vec3 pos = r.origin();
    const Vec3& dir = r.direction();
    int domain_one = -1;
    int domain_two = -1;
    ImageBasedVolumeEngine::Hit hit;
    while (travelled <= max_depth)
    {
        // Evaluate both implicit functions.
        double dist_one = l_field_one->getValue(pos, domain_one);
        double dist_two = l_field_two->getValue(pos, domain_two);
        double dist = fabs(fmax(dist_one, dist_two));
        if(dist<eps)
        {
            // Fill in hit information and return.
            hit.found = true;
			hit.pos = pos;
			hit.normal = (dist_one < dist_two)? l_field_one->getGradient(pos, domain_one) : l_field_two->getGradient(pos, domain_two);
			hit.distance = travelled;
            hit.surface_id = (dist_one < dist_two)? 1 : 0;
            hit.domain = (dist_one < dist_two)? domain_one : domain_two;
            return hit;
        }
        travelled += dist;
        pos += dir * dist;
    }
    return hit;
}

void ImageBasedVolumeEngine::doUpdate()
{
    if ( l_field_one.empty() || l_field_two.empty() )
    {
        return;
    }

    ///// THEN tell everthing is (will be) up to date now
    /// @warning This must be done AFTER updating all inputs
    /// can be done before or after setting up the outputs
    cleanDirty();

    // Initialize accessors.
    // Inputs.
    sofa::defaulttype::BoundingBox bbox_one(l_field_one->l_dofs->f_bbox.getValue());
    auto dof_one = getReadAccessor(*l_field_one->l_dofs->read(sofa::core::VecCoordId::position()));
    sofa::defaulttype::BoundingBox bbox_two(l_field_two->l_dofs->f_bbox.getValue());
    auto dof_two = getReadAccessor(*l_field_two->l_dofs->read(sofa::core::VecCoordId::position()));
    auto res = getReadAccessor(d_resolution);
    double eps = getReadAccessor(d_epsilon);
    // Outputs.
    auto intersections = getWriteAccessor(d_intersections);
    double volume; // = getWriteAccessor(d_volume);
    auto volume_gradients_one = getWriteAccessor(d_volume_gradients_one);
    auto volume_gradients_two = getWriteAccessor(d_volume_gradients_two);

    // Clear outputs.
    intersections.clear();
    volume = 0.0;
    volume_gradients_one.clear();
    volume_gradients_one.resize(dof_one.size());
    volume_gradients_two.clear();
    volume_gradients_two.resize(dof_two.size());

    // Initialize containers.
    Vec3 bbox_bottom, bbox_size, temp;
    sofa::defaulttype::Vec4d barycentric_coordinates;
    sofa::core::topology::BaseMeshTopology::Tetrahedron tetra;

    // Broad phase 
    if (bbox_one.intersect(bbox_two))
    {
        // Construct the AABB of the (potential) interpenetration volume.
        auto bbox = bbox_one.getIntersection(bbox_two);
        bbox_bottom = bbox.minBBox();
        bbox_size = bbox.maxBBox()-bbox_bottom;
    }
    else
    {
        msg_warning() << "The broad phase determined that there is no interpenetration.";
        return;
    }
    // Narrow phase 
    
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
        double max_depth = bbox_size[0]*viewing_direction[0] + bbox_size[1]*viewing_direction[1] + bbox_size[2]*viewing_direction[2]; // pointwise vector multiplication
        // Begin ray casting.
        Vec3 current_line = bbox_bottom;
        for (int i=0; i<res->x(); i++)
        {
            Vec3 current_column = current_line;
            for (int j=0; j<res->y(); j++)
            {
                // Launch ray.
                Ray ray {current_column, viewing_direction};
                ImageBasedVolumeEngine::Hit hit = sphereTracing(ray, eps, max_depth);
                if (hit.found) 
                {
                    // Store interesection point (for display purposes).
                    intersections.push_back(hit.pos);
                    // Accumulate Volume.
                    volume -= pixel_area * hit.distance;
                    // Accumulate volume gradients.
                    if (hit.surface_id)
                    {
                        barycentric_coordinates = l_field_one->getBarycentricCoordinates(hit.pos, hit.domain, dof_one);
                        tetra = l_field_one->l_topology->getTetrahedron(hit.domain );
                        for (unsigned int k=0; k<4; k++)
                        {
                            // TODO: ask Damien about pointwise vector multiplication
                            temp[0] = dof_one[tetra[k]][0] * viewing_direction[0];
                            temp[1] = dof_one[tetra[k]][1] * viewing_direction[1];
                            temp[2] = dof_one[tetra[k]][2] * viewing_direction[2];
                            volume_gradients_one[hit.domain] += -1 * pixel_area * barycentric_coordinates[k] * temp;
                        }
                    }
                    else
                    {
                        barycentric_coordinates = l_field_two->getBarycentricCoordinates(hit.pos, hit.domain, dof_two);
                        tetra = l_field_two->l_topology->getTetrahedron(hit.domain);
                        for (unsigned int k=0; k<4; k++)
                        {
                            // TODO: ask Damien about pointwise vector multiplication
                            temp[0] = dof_two[tetra[k]][0] * viewing_direction[0];
                            temp[1] = dof_two[tetra[k]][1] * viewing_direction[1];
                            temp[2] = dof_two[tetra[k]][2] * viewing_direction[2];
                            volume_gradients_two[hit.domain] += -1 * pixel_area * barycentric_coordinates[k] * temp;
                        }
                    }                    
                    // Launch ray from the oposite side of the bounding box.
                    ray.setOrigin(current_column+max_depth*viewing_direction);
                    ray.setDirection(viewing_direction * -1);
                    ImageBasedVolumeEngine::Hit hit = sphereTracing(ray, eps, max_depth);
                    if (hit.found)
                    {
                        // Store interesection point (for display purposes).
                        intersections.push_back(hit.pos);
                        // Accumulate Volume.
                        volume += pixel_area * (max_depth-hit.distance);
                        // Accumulate volume gradients.
                        if (hit.surface_id)
                        {
                            barycentric_coordinates = l_field_one->getBarycentricCoordinates(hit.pos, hit.domain, dof_one);
                            tetra = l_field_one->l_topology->getTetrahedron(hit.domain);
                            for (unsigned int k=0; k<4; k++)
                            {
                                // TODO: ask Damien about pointwise vector multiplication
                                temp[0] = dof_one[tetra[k]][0] * viewing_direction[0];
                                temp[1] = dof_one[tetra[k]][1] * viewing_direction[1];
                                temp[2] = dof_one[tetra[k]][2] * viewing_direction[2];
                                volume_gradients_one[hit.domain] += pixel_area * barycentric_coordinates[k] * temp;
                            }
                        }
                        else
                        {
                            barycentric_coordinates = l_field_two->getBarycentricCoordinates(hit.pos, hit.domain, dof_two);
                            tetra = l_field_two->l_topology->getTetrahedron(hit.domain);
                            for (unsigned int k=0; k<4; k++)
                            {
                                // TODO: ask Damien about pointwise vector multiplication
                                temp[0] = dof_two[tetra[k]][0] * viewing_direction[0];
                                temp[1] = dof_two[tetra[k]][1] * viewing_direction[1];
                                temp[2] = dof_two[tetra[k]][2] * viewing_direction[2];
                                volume_gradients_two[hit.domain] += pixel_area * barycentric_coordinates[k] * temp;
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

    d_volume.beginEdit();
    d_volume.setValue(volume);
    d_volume.endEdit();

    return;
}

void ImageBasedVolumeEngine::draw(const sofa::core::visual::VisualParams* params)
{
    doUpdate(); // CHEAT
}

} /// namespace sofa::component::engine