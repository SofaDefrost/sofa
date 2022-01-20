#include <sofa/core/ObjectFactory.h>
#include "ImageBasedVolumeEngineOgl.h"

/*

AFFICHER NOMBRE D'ITERATION 

COMPARER AVEC RAY MARCHING SUR CPU

AFFICHER SI ON A TOUCHER UN TETRAHEDRE OU NON

IMPLEMENTER EN CPU L'ALGO EN UTILISANT L'INTERSECTION RAYON-TETRAHEDRE

AFFICHER LES RESULTATS UNE ITERATION A LA FOIS
AFFICHER LES RESULTATS UNE ITERATION A LA FOIS
AFFICHER LES RESULTATS UNE ITERATION A LA FOIS
AFFICHER LES RESULTATS UNE ITERATION A LA FOIS

*/

namespace sofa::component::engine
{
using sofa::defaulttype::Ray;
using sofa::helper::getReadAccessor;
using sofa::helper::getWriteAccessor;
using sofaimplicitfield::DisplacementField;
using sofa::defaulttype::Matrix4;

/// Register in the Factory.
static int ImageBasedVolumeEngineOglClass = core::RegisterObject("Mono-volume setting. GPU implementation using OpenGL.").add< ImageBasedVolumeEngineOgl >();

// Inner class.
class ImageBasedVolumeEngineOgl::InternalData
{
public:
    // Attached shader:
    sofa::component::visualmodel::OglShader* shader;
    // Frame buffer:
    GLuint outFrameBufferObj, outDepthStencilBuffer;
    unsigned int width;
    unsigned int height;
    // Textures:
    GLuint dofsBuffer, dofsTexture;
    std::vector<sofa::type::Vec3f> dofs;
    GLuint restBuffer, restTexture;
    std::vector<sofa::type::Vec3f> rest;
    // Input buffers:
    unsigned int nbPassedIndices;
    GLuint in_tetrahedraIndicesBuffer;
    std::vector<int> in_tetrahedraIndices;
    // Output buffers:
    GLenum buffs[5] = {GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1, GL_COLOR_ATTACHMENT2, GL_COLOR_ATTACHMENT3, GL_COLOR_ATTACHMENT4};
    GLuint hitDataDepthBuffer, hitDataBarycentricWeightsBuffer, hitDataTetrahedronIDBuffer, hitDataNormalsBuffer, hitDataPosBuffer;
    std::vector<float> front_hitDataDepth;
    std::vector<sofa::type::Vec4f> front_hitDataBarycentricWeights;
    std::vector<int> front_hitDataTetrahedronID;
    std::vector<sofa::type::Vec3f> front_hitDataNormals;
    std::vector<sofa::type::Vec3f> front_hitDataPos;
    std::vector<float> back_hitDataDepth;
    std::vector<sofa::type::Vec4f> back_hitDataBarycentricWeights;
    std::vector<int> back_hitDataTetrahedronID;
    std::vector<sofa::type::Vec3f> back_hitDataNormals;
    std::vector<sofa::type::Vec3f> back_hitDataPos;
    // Storage:
    float viewingDirectionArray[3];
    // Boolean status trackers:
    bool isInited = false;
    bool isStaticDataSet = false;

    void init(sofa::component::visualmodel::OglShader* targetShader, unsigned int width, unsigned int height)
    {
        std::cout << "ImageBasedVolumeEngineOgl::InternalData init()" << std::endl;
        // CPU-cide initialization:
        this->width = width;
        this->height = height;
        // Resize output buffers.
        front_hitDataDepth.resize(width * height);
        back_hitDataDepth.resize(width * height);
        front_hitDataBarycentricWeights.resize(width * height);
        back_hitDataBarycentricWeights.resize(width * height);
        front_hitDataTetrahedronID.resize(width * height);
        back_hitDataTetrahedronID.resize(width * height);
        front_hitDataNormals.resize(width * height);
        back_hitDataNormals.resize(width * height);
        front_hitDataPos.resize(width * height);
        back_hitDataPos.resize(width * height);
        
        // GPU-side initialization.
        shader = targetShader;
        // Initialize frame buffer.
        glGenFramebuffers(1, &outFrameBufferObj);
        glBindFramebuffer(GL_FRAMEBUFFER, outFrameBufferObj);
        // Initialize the frame buffer's depth and stencil buffer:
        glGenRenderbuffers(1, &outDepthStencilBuffer);
        glBindRenderbuffer(GL_RENDERBUFFER, outDepthStencilBuffer);
        glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH24_STENCIL8, width, height);
        glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_RENDERBUFFER, outDepthStencilBuffer);
        glBindRenderbuffer(GL_RENDERBUFFER, 0);
        // Initialize the frame buffer's outputs:
        // 1) Ray-marching hit's depth.
        glGenTextures(1, &hitDataBarycentricWeightsBuffer);
        glBindTexture(GL_TEXTURE_2D, hitDataBarycentricWeightsBuffer);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, width, height, 0, GL_RED, GL_FLOAT, NULL);
        glFramebufferTexture2D(GL_FRAMEBUFFER, buffs[0], GL_TEXTURE_2D, hitDataBarycentricWeightsBuffer, 0);
        // 2) Ray-marching hit's barycentric coordinates.
        glGenTextures(1, &hitDataBarycentricWeightsBuffer); 
        glBindTexture(GL_TEXTURE_2D, hitDataBarycentricWeightsBuffer);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, width, height, 0, GL_RGBA, GL_FLOAT, NULL);
        glFramebufferTexture2D(GL_FRAMEBUFFER, buffs[1],  GL_TEXTURE_2D, hitDataBarycentricWeightsBuffer, 0); 
        // 3) Ray-marching hit's domain.
        glGenTextures(1, &hitDataTetrahedronIDBuffer); 
        glBindTexture(GL_TEXTURE_2D, hitDataTetrahedronIDBuffer);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_R8I, width, height, 0, GL_RED_INTEGER, GL_INT, NULL);
        glFramebufferTexture2D(GL_FRAMEBUFFER, buffs[2],  GL_TEXTURE_2D, hitDataTetrahedronIDBuffer, 0); 
        // 4) Ray-marching hit's normal.
        glGenTextures(1, &hitDataNormalsBuffer);
        glBindTexture(GL_TEXTURE_2D, hitDataNormalsBuffer);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB32F, width, height, 0, GL_RGB, GL_FLOAT, NULL);
        glFramebufferTexture2D(GL_FRAMEBUFFER, buffs[3],  GL_TEXTURE_2D, hitDataNormalsBuffer, 0);
        // 5) Ray-marching hit's position.
        glGenTextures(1, &hitDataPosBuffer);
        glBindTexture(GL_TEXTURE_2D, hitDataPosBuffer);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB32F, width, height, 0, GL_RGB, GL_FLOAT, NULL);
        glFramebufferTexture2D(GL_FRAMEBUFFER, buffs[4],  GL_TEXTURE_2D, hitDataPosBuffer, 0);
        glBindTexture(GL_TEXTURE_2D, 0);

        // Define the frame buffer's outputs.
        glDrawBuffers(5, buffs);
        // Disable clamp outputs to [0,1].
        glClampColor(GL_CLAMP_READ_COLOR, GL_FALSE);
        // Clear.
        glClearDepthf(1.0f); // needs to be 1 because we use GL_LESS!
        glClearColor(-1.0f, -1.0f, -1.0f, -1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

        // Check frame buffer sanity:
        if(glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
        {
            std::cout << "Failed frame buffer sanity check" << std::endl;
            return;
        }
        glBindFramebuffer(GL_FRAMEBUFFER, 0);

        // Initialize input vertex buffers.
        glGenBuffers(1, &in_tetrahedraIndicesBuffer);
        // Initialize uniform texture buffers.
        // ... current positons:
        glGenBuffers(1, &dofsBuffer);
        glGenTextures(1, &dofsTexture);
        // ... rest positions:     
        glGenBuffers(1, &restBuffer);
        glGenTextures(1, &restTexture);

        isInited = true;
        return;
    }

    void setStaticData(sofa::helper::ReadAccessor<sofa::helper::vector<sofa::defaulttype::Vec3d>>& restPositions, 
                       const sofa::core::topology::BaseMeshTopology::SeqTetrahedra& tetrahedra,
                       /*std::vector<sofa::core::topology::BaseMeshTopology::Tetrahedron>& tetrahedra,*/
                       float smallStep, float epsilon, float h)
    {
        std::cout << "ImageBasedVolumeEngineOgl::InternalData setStaticData()" << std::endl;
        // Fill rest positions
        rest.resize(restPositions.size());
        for(sofa::Index i=0; i<restPositions.size(); i++)
        {
            rest[i] = (sofa::type::Vec3f) restPositions[i];
        }
        // Fill tetrahedra indices
        for(auto tetra : tetrahedra)
        {
            in_tetrahedraIndices.push_back(tetra[0]);
            in_tetrahedraIndices.push_back(tetra[1]);
            in_tetrahedraIndices.push_back(tetra[2]);
            in_tetrahedraIndices.push_back(tetra[3]);
        }
        nbPassedIndices = in_tetrahedraIndices.size();

        // Access uniform data fields:
        GLuint u_epsilon = shader->getUniform(shader->getCurrentIndex(), "epsilon");
        GLuint u_h = shader->getUniform(shader->getCurrentIndex(), "h");
        GLuint u_smallStep = shader->getUniform(shader->getCurrentIndex(), "smallStep");

        // Start attached shader.
        shader->start(); //glBindFramebuffer(GL_FRAMEBUFFER, outFrameBufferObj);
        // Bind and fill rest positions.
        glBindBuffer(GL_TEXTURE_BUFFER, restBuffer);
        glBufferData(GL_TEXTURE_BUFFER, sizeof(sofa::type::Vec3f)*rest.size(), rest.data(), GL_STATIC_DRAW);
        glBindBuffer(GL_TEXTURE_BUFFER, 0);

        // Bind and fill indices.
        glBindBuffer(GL_ARRAY_BUFFER, in_tetrahedraIndicesBuffer);
        glBufferData(GL_ARRAY_BUFFER, sizeof(int) * in_tetrahedraIndices.size(), in_tetrahedraIndices.data(), GL_STATIC_DRAW);
        glBindBuffer(GL_ARRAY_BUFFER, 0);

        // Pass uniform data.
        glUniform1f(u_epsilon, epsilon);
        glUniform1f(u_h, h);
        glUniform1f(u_smallStep, smallStep);

        // Stop attached shader.
        shader->stop();//glBindFramebuffer(GL_FRAMEBUFFER, 0);

        isStaticDataSet = true;
        return;
    }

    void updateCurrentPositionTexture(sofa::helper::ReadAccessor<sofa::helper::vector<sofa::defaulttype::Vec3d>>& currentPositions)
    {
        std::cout << "ImageBasedVolumeEngineOgl::InternalData updateCurrentPositionTexture()" << std::endl;
        // Fill dofs positions.
        dofs.resize(currentPositions.size());
        for(sofa::Index i=0; i<currentPositions.size(); i++)
        {
            dofs[i] = (sofa::type::Vec3f) currentPositions[i];
        }

        // Start attached shader.
        shader->start(); //glBindFramebuffer(GL_FRAMEBUFFER, outFrameBufferObj);
        // Bind and fill current positions.
        glBindBuffer(GL_TEXTURE_BUFFER, dofsBuffer);
        glBufferData(GL_TEXTURE_BUFFER, sizeof(sofa::type::Vec3f)*dofs.size(), dofs.data(), GL_DYNAMIC_DRAW);
        glBindBuffer(GL_TEXTURE_BUFFER, 0);
        // Stop attached shader.
        shader->stop(); //glBindFramebuffer(GL_FRAMEBUFFER, 0);
        return;
    }

    void updateUniformDataFields(float* modelMatArray, float* viewMatArray, float* projMatArray, Vec3 viewingDirection, float zNear, float zFar)
    {
        std::cout << "ImageBasedVolumeEngineOgl::InternalData updateUniformDataFields()" << std::endl;

        // Type conversions.
        viewingDirectionArray[0] = float(viewingDirection[0]);
        viewingDirectionArray[1] = float(viewingDirection[1]);
        viewingDirectionArray[2] = float(viewingDirection[2]);
        // Set values.
        shader->setMatrix4(shader->getCurrentIndex(), "modelMatrix", 1, false, modelMatArray);
        shader->setMatrix4(shader->getCurrentIndex(), "viewMatrix", 1, false, viewMatArray);
        shader->setMatrix4(shader->getCurrentIndex(), "projectionMatrix", 1, false, projMatArray);
        //shader->setFloatVector3(shader->getCurrentIndex(), "viewingDirection", 1, viewingDirectionArray);
        shader->setFloat(shader->getCurrentIndex(), "zNear", zNear);
        shader->setFloat(shader->getCurrentIndex(), "zFar", zFar);
        return;
    }

    void renderAndRead()
    {
        std::cout << "ImageBasedVolumeEngineOgl::InternalData renderAndRead()" << std::endl;
        // Attach FBO.
        glBindFramebuffer(GL_FRAMEBUFFER, outFrameBufferObj);

        //
        glClearColor(-1.0f, -1.0f, -1.0f, -1.0f);
        glClearDepthf(1.0f);
        glClearStencil(0);

        // Enable early testings.
        glEnable(GL_DEPTH_TEST);
        glEnable(GL_STENCIL_TEST);

        // Enable rendering of 'back-facing' polygons.
        glDisable(GL_CULL_FACE);

        // Define viewport.
        glViewport(0, 0, width, height);

        // Uniform addresses.
        GLuint u_dofSampler = shader->getUniform(shader->getCurrentIndex(), "dofBuffer");
        GLuint u_restSampler = shader->getUniform(shader->getCurrentIndex(), "restBuffer");
        
        GLuint u_viewingDirection = shader->getUniform(shader->getCurrentIndex(), "viewingDirection");

        // Enable/bind attributes.
        glEnableVertexAttribArray(0);
        glBindBuffer(GL_ARRAY_BUFFER, in_tetrahedraIndicesBuffer);
        glVertexAttribIPointer(0, 1, GL_INT, 0, 0); // layout(location=0) in int

        // Attach shader.
        shader->start();

        // Activate/Bind uniform samplers.
        // ... current positions
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_BUFFER, dofsTexture);
        glTexBuffer(GL_TEXTURE_BUFFER, GL_RGB32F, dofsBuffer);
        glUniform1i(u_dofSampler, 0);
        // ... rest positions
        glActiveTexture(GL_TEXTURE1);
        glBindTexture(GL_TEXTURE_BUFFER, restTexture);
        glTexBuffer(GL_TEXTURE_BUFFER, GL_RGB32F, restBuffer);
        glUniform1i(u_restSampler, 1);

        // CLEAR EVERYTHING!
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
        
        // Draw front to back.
        glDepthFunc(GL_LESS);

        glStencilFunc(GL_ALWAYS, 1, 0xFF); // Set any stencil to 1
        glStencilOp(GL_KEEP, GL_KEEP, GL_REPLACE);
        glStencilMask(0xFF); // Write to stencil buffer

        // Set ray-direction.
        glUniform3f(u_viewingDirection, viewingDirectionArray[0], viewingDirectionArray[1], viewingDirectionArray[2]);

        glDrawArrays(GL_LINES_ADJACENCY, 0, nbPassedIndices); 
        // "GL_LINES_ADJACENCY" feels like cheating. But it allows us to handle 4 vertices at a time in the geometry shader.

        // Read results:
        // ..ray-marching hit's depth.
        glReadBuffer(buffs[0]);
        glReadPixels( 0, 0, width, height, GL_RED,  GL_FLOAT, front_hitDataDepth.data());
        // ... ray-marching hit's barycentric coordinates.
        glReadBuffer(buffs[1]);
        glReadPixels( 0, 0, width, height, GL_RGBA,  GL_FLOAT, front_hitDataBarycentricWeights.data());
        // ... ray-marching hit's domain.
        glReadBuffer(buffs[2]);
        glReadPixels( 0, 0, width, height, GL_RED_INTEGER,  GL_INT, front_hitDataTetrahedronID.data());
        // ... ray-marching hit's normal.
        glReadBuffer(buffs[3]);
        glReadPixels( 0, 0, width, height, GL_RGB,  GL_FLOAT, front_hitDataNormals.data());
        // ... ray-marching hit's position.
        glReadBuffer(buffs[4]);
        glReadPixels( 0, 0, width, height, GL_RGB,  GL_FLOAT, front_hitDataPos.data());

        // Draw back to front.
        glDepthFunc(GL_GEQUAL);  

        glStencilFunc(GL_EQUAL, 1, 0xFF); // Only consider fragments already belonging to the pixel-coverage-set
        glStencilMask(0x00); // Don't write anything to stencil buffer

        // Set ray-direction.
        glUniform3f(u_viewingDirection, -1*viewingDirectionArray[0], -1*viewingDirectionArray[1], -1*viewingDirectionArray[2]);

        glDrawArrays(GL_LINES_ADJACENCY, 0, nbPassedIndices); 
        // "GL_LINES_ADJACENCY" feels like cheating. But it allows us to handle 4 vertices at a time in the geometry shader.

        // Read results:
        // ... ray-marching hit's depth.
        glReadBuffer(buffs[0]);
        glReadPixels( 0, 0, width, height, GL_RED,  GL_FLOAT, back_hitDataDepth.data());
        // ... ray-marching hit's barycentric coordinates.
        glReadBuffer(buffs[1]);
        glReadPixels( 0, 0, width, height, GL_RGBA,  GL_FLOAT, back_hitDataBarycentricWeights.data());
        // ... ray-marching hit's domain.
        glReadBuffer(buffs[2]);
        glReadPixels( 0, 0, width, height, GL_RED_INTEGER,  GL_INT, back_hitDataTetrahedronID.data());
        // ... ray-marching hit's normal.
        glReadBuffer(buffs[3]);
        glReadPixels( 0, 0, width, height, GL_RGB,  GL_FLOAT, back_hitDataNormals.data());
        // ... ray-marching hit's position.
        glReadBuffer(buffs[4]);
        glReadPixels( 0, 0, width, height, GL_RGB,  GL_FLOAT, back_hitDataPos.data());

        // CLEAR EVERYTHING!
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

        // Deactivate/Unbind uniform samplers.
        glActiveTexture(GL_TEXTURE1);
        glBindTexture(GL_TEXTURE_BUFFER, 0);
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_BUFFER, 0);

        // Disable/unbind attributes.
        glDisableVertexAttribArray(0);
        glBindBuffer(GL_ARRAY_BUFFER, 0);

        // Disable early testings.
        glDisable(GL_DEPTH_TEST);
        glDisable(GL_STENCIL_TEST);

        // Disable rendering of 'back-facing' polygons.
        //glEnable(GL_CULL_FACE);

        // Detach shader and FBO.
        shader->stop();
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        return;
    }
};

// Constructor.
ImageBasedVolumeEngineOgl::ImageBasedVolumeEngineOgl():
    // Inputs:
    l_field_one(initLink("field_one", "The first colliding scalar field.")),
    l_shader_one(initLink("shader_one", "The shader to use for sampling the first scallar field.")),
    l_field_two(initLink("field_two", "The second colliding scalar field.")),
    l_shader_two(initLink("shader_two", "The shader to use for sampling the first scallar field.")),
    d_resolution(initData(&d_resolution, Vec2i{20,20}, "resolution", "The amount of samples per visual axis.")),
    d_epsilon(initData(&d_epsilon, 0.01, "epsilon", "The tolerance allowed when evaluating an implicit surface.")),
    // Outputs:
    d_intersections(initData(&d_intersections, "intersections", "The intersection points' locations.")),
    d_volume(initData(&d_volume, "volume", "The evaluated interpenetration volume.")),
    d_volume_gradients_one(initData(&d_volume_gradients_one, "volume_gradients_one", "The evaluated interpenetration volume' gradients w.r.t. the DOFs of the first scallar field.")),
    d_volume_gradients_two(initData(&d_volume_gradients_two, "volume_gradients_two", "The evaluated interpenetration volume' gradients w.r.t. the DOFs of the second scallar field."))
{
    // Attach internal class instance.
    data_one.reset(new InternalData());
    data_two.reset(new InternalData());
    // Declare outputs.
    addOutput(&d_intersections);
    addOutput(&d_volume);
    addOutput(&d_volume_gradients_one);
    addOutput(&d_volume_gradients_two);
}

// Public.
void ImageBasedVolumeEngineOgl::init()
{
    setDirtyValue();
}

void ImageBasedVolumeEngineOgl::reinit()
{
    setDirtyValue();
    update();
}

void ImageBasedVolumeEngineOgl::doUpdate()
{
    // Abbort early calls.
    if ( l_field_one.empty() || l_field_two.empty() )
    {
        return;
    }
    
    if( !l_shader_one.get() || !l_shader_two.get() )
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
    double volume = getWriteAccessor(d_volume);
    auto volume_gradients_one = getWriteAccessor(d_volume_gradients_one);
    auto volume_gradients_two = getWriteAccessor(d_volume_gradients_two);

    // Clear outputs.
    intersections.clear();
    volume = 0.0;
    volume_gradients_one.clear();
    volume_gradients_one.resize(dof_one.size());
    volume_gradients_two.clear();
    volume_gradients_two.resize(dof_two.size());

    /** 
     * Broad phase collision detection
    **/
    Vec3 bboxInter_bot, bboxInter_top, bboxInter_size;
    Vec3 bboxUnion_bot, bboxUnion_top, bboxUnion_size;
    // Test AABB intersection.
    if (bbox_one.intersect(bbox_two))
    {
        // Construct the AABBs of the (potential) interpenetration volume.
        // ... intersection of the colliding bboxes:
        auto bbox = bbox_one.getIntersection(bbox_two);
        bboxInter_bot = bbox.minBBox();
        bboxInter_top = bbox.maxBBox();
        bboxInter_size = bboxInter_top-bboxInter_bot;
        // ... "" union "" of the colliding bboxes:
        bboxUnion_bot = {fmin(bbox_one.minBBox()[0], bbox_two.minBBox()[0]), fmin(bbox_one.minBBox()[1], bbox_two.minBBox()[1]), fmin(bbox_one.minBBox()[2], bbox_two.minBBox()[2])};
        bboxUnion_top = {fmax(bbox_one.maxBBox()[0], bbox_two.maxBBox()[0]), fmax(bbox_one.maxBBox()[1], bbox_two.maxBBox()[1]), fmax(bbox_one.maxBBox()[2], bbox_two.maxBBox()[2])};
        bboxUnion_size = bboxUnion_top-bboxUnion_bot;
    }
    else
    {
        msg_warning() << "The broad phase determined that there is no interpenetration.";
        return;
    }

    msg_warning() << "AABB n: " << bboxInter_bot << " -> " << bboxInter_top;
    msg_warning() << "size: " << bboxInter_size;
    msg_warning() << "AABB u: " << bboxUnion_bot << " -> " << bboxUnion_top;
    msg_warning() << "size: " << bboxUnion_size;

    /** 
     * Narrow phase collision detection
    **/
    Vec3 tempVec3;
    sofa::core::topology::BaseMeshTopology::Tetrahedron tetra;

    // Initialize internal data.
    if(!data_one->isInited)
    {
        data_one->init(l_shader_one, res->x(), res->y());
    }
    if(!data_two->isInited)
    {
        data_two->init(l_shader_two, res->x(), res->y());
    }
    
    // Set-up internal data's 'static' buffers.
    if(!data_one->isStaticDataSet)
    {
        auto restDofs = getReadAccessor(*l_field_one->l_dofs->read(sofa::core::VecCoordId::restPosition()));
        data_one->setStaticData(restDofs, l_field_one->l_topology->getTetrahedra(), 0.001f, float(eps), 0.00001f); // I'd want to pass float(l_field_one->d_epsilon.getValue() but it's ~protected~
    }
    if(!data_two->isStaticDataSet)
    {
        auto restDofs = getReadAccessor(*l_field_two->l_dofs->read(sofa::core::VecCoordId::restPosition()));
        data_two->setStaticData(restDofs, l_field_two->l_topology->getTetrahedra(), 0.001f, float(eps), 0.00001f); // I'd want to pass float(l_field_two->d_epsilon.getValue() but it's ~protected~
    }

    // Pass current DOFs positions to internal data.
    data_one->updateCurrentPositionTexture(dof_one);
    data_two->updateCurrentPositionTexture(dof_two);

    // Iterate over the faces of the AABB, pairing the opposite faces together.
    std::vector<sofa::defaulttype::Vec2i> planes {{1,2}, {0,1}, {0,2}};

    for (unsigned int plane_it=0; plane_it<3; plane_it++)
    {
        // Compute area of each pixel.
        double pixel_area = (bboxInter_size[planes[plane_it][0]] / res->x()) * (bboxInter_size[planes[plane_it][1]] / res->y());

        // Define "camera" basis.
        Vec3 viewingDirection {1.0, 1.0, 1.0};
        viewingDirection[planes[plane_it][0]] = 0.0;
        viewingDirection[planes[plane_it][1]] = 0.0;
        Vec3 widthDirection {0.0, 0.0, 0.0};
        widthDirection[planes[plane_it][0]] = 1.0;
        Vec3 heightDirection {0.0, 0.0, 0.0};
        heightDirection[planes[plane_it][1]] = 1.0;
        
        // Define model matrix
        float modelMatArray[16] = {0.0};
        modelMatArray[0] = widthDirection[0];
        modelMatArray[4] = widthDirection[1];
        modelMatArray[8] = widthDirection[2];
        modelMatArray[1] = heightDirection[0];
        modelMatArray[5] = heightDirection[1];
        modelMatArray[9] = heightDirection[2];
        modelMatArray[2] = viewingDirection[0];
        modelMatArray[6] = viewingDirection[1];
        modelMatArray[10] = viewingDirection[2];
        modelMatArray[15] = 1;
    
        // Define view matrix.
        float viewMatArray[16] = {0.0};
        viewMatArray[0] = 1;
        viewMatArray[5] = 1;
        viewMatArray[10] = 1;
        viewMatArray[15] = 1;
        viewMatArray[12] = plane_it == 0 ? -1*(bboxUnion_bot[0]) : -1*(bboxInter_bot[0]);
        viewMatArray[13] = plane_it == 1 ? -1*(bboxUnion_bot[1]) : -1*(bboxInter_bot[1]);
        viewMatArray[14] = plane_it == 2 ? -1*(bboxUnion_bot[2]) : -1*(bboxInter_bot[2]);

        // Define orthogonal projection matrix.
        // 2/(right-left)   0               0               -(right+left) / (right-left)
        // 0                2/(top-bottom)  0               -(top+bottom) / (top-bottom)
        // 0                0               -2/(far-near)    -(far+near) / (far-near)
        // 0                0               0               1
        float projectionMatArray[16] = {0.0};
        projectionMatArray[0] = 2 / bboxInter_size[planes[plane_it][0]];
        projectionMatArray[5] = 2 / bboxInter_size[planes[plane_it][1]];
        projectionMatArray[10] = 2 / bboxUnion_size[plane_it];
        projectionMatArray[12] = -1;
        projectionMatArray[13] = -1;
        projectionMatArray[14] = -1;
        projectionMatArray[15] = 1;

        // Update internal data's uniform fields.
        data_one->updateUniformDataFields(modelMatArray, viewMatArray, projectionMatArray, viewingDirection,
                                          float(bboxUnion_bot[plane_it]), float(bboxUnion_top[plane_it]));
        data_two->updateUniformDataFields(modelMatArray, viewMatArray, projectionMatArray, viewingDirection,
                                          float(bboxUnion_bot[plane_it]), float(bboxUnion_top[plane_it]));

        // Sample intersection Volume.
        data_one->renderAndRead();
        data_two->renderAndRead();

        // Accumulate results.
        for (int i=0; i<res->x()*res->y(); i++)
        {
            // (front) Check that the same pixel covers both objects.
            if (data_one->front_hitDataDepth[i] != -1 && data_two->front_hitDataDepth[i] != -1)
            {
                // Identify to which surface the sample of the interpenetration volume bolongs to.
                // ... (front)
                if (data_one->front_hitDataDepth[i] > data_two->front_hitDataDepth[i])
                {
                    // Fill d_intersections.
                    intersections.push_back(data_one->front_hitDataPos[i]);
                    // V.
                    volume -= pixel_area * double(data_one->front_hitDataDepth[i]);
                    // dV.
                    tetra = l_field_one->l_topology->getTetrahedron(data_one->front_hitDataTetrahedronID[i]);
                    for (unsigned int j=0; j<4; j++)
                    {
                        tempVec3[0] = dof_one[tetra[j]][0] * viewingDirection[0];
                        tempVec3[1] = dof_one[tetra[j]][1] * viewingDirection[1];
                        tempVec3[2] = dof_one[tetra[j]][2] * viewingDirection[2];
                        volume_gradients_one[data_one->front_hitDataTetrahedronID[i]] += (-1) * pixel_area * data_one->front_hitDataBarycentricWeights[i][j] * tempVec3;
                    }
                }
                else
                {
                    // Fill d_intersections.
                    intersections.push_back(data_two->front_hitDataPos[i]);
                    // V.
                    volume -= pixel_area * double(data_two->front_hitDataDepth[i]);
                    // dV.
                    tetra = l_field_two->l_topology->getTetrahedron(data_two->front_hitDataTetrahedronID[i]);
                    for (unsigned int j=0; j<4; j++)
                    {
                        tempVec3[0] = dof_two[tetra[j]][0] * viewingDirection[0];
                        tempVec3[1] = dof_two[tetra[j]][1] * viewingDirection[1];
                        tempVec3[2] = dof_two[tetra[j]][2] * viewingDirection[2];
                        volume_gradients_two[data_two->front_hitDataTetrahedronID[i]] += (-1) * pixel_area * data_two->front_hitDataBarycentricWeights[i][j] * tempVec3;
                    }
                }
                // ... (back)
                if (data_one->back_hitDataDepth[i] < data_two->back_hitDataDepth[i])
                {
                    // fill d_intersections !
                    intersections.push_back(data_one->back_hitDataPos[i]);
                    // V
                    volume += pixel_area * double(data_one->back_hitDataDepth[i]);
                    // dV
                    tetra = l_field_one->l_topology->getTetrahedron(data_one->back_hitDataTetrahedronID[i]);
                    for (unsigned int j=0; j<4; j++)
                    {
                        tempVec3[0] = dof_one[tetra[j]][0] * viewingDirection[0];
                        tempVec3[1] = dof_one[tetra[j]][1] * viewingDirection[1];
                        tempVec3[2] = dof_one[tetra[j]][2] * viewingDirection[2];
                        volume_gradients_one[data_one->back_hitDataTetrahedronID[i]] += pixel_area * data_one->back_hitDataBarycentricWeights[i][j] * tempVec3;
                    }
                }
                else
                {

                    std::cout << 'd' << std::endl;

                    // fill d_intersections !
                    intersections.push_back(data_two->back_hitDataPos[i]);
                    // V
                    volume += pixel_area * double(data_two->back_hitDataDepth[i]);
                    // dV
                    tetra = l_field_two->l_topology->getTetrahedron(data_two->back_hitDataTetrahedronID[i]);

                    std::cout << "d bis" << std::endl;

                    std::cout << "data_two->back_hitDataTetrahedronID[i]: " << data_two->back_hitDataTetrahedronID[i] << std::endl;

                    for (unsigned int j=0; j<4; j++)
                    {
                        tempVec3[0] = dof_two[tetra[j]][0] * viewingDirection[0];
                        tempVec3[1] = dof_two[tetra[j]][1] * viewingDirection[1];
                        tempVec3[2] = dof_two[tetra[j]][2] * viewingDirection[2];
                        volume_gradients_two[data_two->back_hitDataTetrahedronID[i]] += pixel_area * data_two->back_hitDataBarycentricWeights[i][j] * tempVec3;
                    }

                    std::cout << "d bis bis" << std::endl;
                }
            }
        }
        msg_warning() << "Volume: " << volume;
    }
    volume /= 3.0;

    msg_warning() << "Final volume: " << volume;

    d_volume.beginEdit();
    d_volume.setValue(volume);
    d_volume.endEdit();

    return;
}

} /// namespace sofa::component::engine