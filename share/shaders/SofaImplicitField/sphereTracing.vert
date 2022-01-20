#version 460
// Vertex shadder.

// Belongs to the OpenGL pipeline:
// sphereTracing_v1.vert -> sphereTracing_v1.geom -> sphereTracing_v1.frag

layout (location = 0) in int vertexIDin;

//out vec4 gl_Position;
out int vertexID; 

uniform mat4 modelMatrix;
uniform mat4 viewMatrix;
uniform mat4 projectionMatrix;
uniform samplerBuffer dofBuffer;

void main()
{
    // Project vertices in the normalized space.
    gl_Position = projectionMatrix * viewMatrix * modelMatrix * vec4(texelFetch(dofBuffer, vertexIDin).xyz, 1);
    // Push through vertices IDs.
    vertexID = vertexIDin;
}