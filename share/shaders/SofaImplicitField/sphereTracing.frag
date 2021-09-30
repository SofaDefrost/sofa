#version 460
// Fragment shadder.

// Belongs to the OpenGL pipeline:
// sphereTracing_v1.vert -> sphereTracing_v1.geom -> sphereTracing_v1.frag

// Before usage this shadder need to be completed with a function float getValue(in vec3 pos){ ... } 
// describing the sampled Signed Distance Function. This can be achieved by writting into the shadder
// as a text stream. Failure to do so will result in a compilation error.

// in in gl_PrimitiveID;
// in bool gl_FrontFacing;
// in vec2 gl_FragCoord;
in vec3 rayOrigin;
flat in vec3 rayDirection;
flat in mat3 restBasisT;
flat in vec3 restBasisAncor;
flat in vec3[4] dofsSummit;

layout(location = 0) out float out_depth;
layout(location = 1) out vec4 out_barycentricWeights;
layout(location = 2) out int out_tetrahedronID;
layout(location = 3) out vec3 out_normal;
layout(location = 4) out vec3 out_pos;

uniform float epsilon;
uniform float h;
uniform float zNear;
uniform float zFar;

float original_depth(float d)
{
    return d * (zFar-zNear) + zNear;
}

float getValue(in vec3 pos)
{
    // This function needs to be defined according to the scene objects.
    // Since shadders are nothing more than text streams, this can be done by 
    // concatenating the computation of the signed distance function in this 
    // dummy function in the initialization of the SOFA scene.
}

vec3 getGradient(in vec3 in_pos)
{
    vec3 normal;
    // Evaluate point.
    float v = getValue(in_pos);

    // Evaluate displaced point:
    in_pos[0] += h;
    normal[0] = getValue(in_pos);
    in_pos[0] -= h;
    in_pos[1] += h;
    normal[1] = getValue(in_pos);
    in_pos[1] -= h;
    in_pos[2] += h;
    normal[2] = getValue(in_pos);
    in_pos[2] -= h;
    
    // Finite difference.
    normal[0] = (normal[0]-v)/h;
    normal[1] = (normal[1]-v)/h;
    normal[2] = (normal[2]-v)/h;
    return normalize(normal);
}


void main()
{
    // If the ray is pointed outward with respect to the normal of the primitive, abbort imediatly.
    if (!gl_FrontFacing)
    {
        discard;
    }

    // Begin sphere tracing.
    vec3 pos = rayOrigin;
    float travelled = 0;
    float dist;
    //vec3 barycentricCoordinates = restBasisT * (pos - restBasisAncor);
    vec3 barycentricCoordinates = {0,0,0};

    // Ray march for as long as the current position is inside the target domain.
    while (barycentricCoordinates.x + barycentricCoordinates.y + barycentricCoordinates.z <= 1 && barycentricCoordinates.x >= 0 && barycentricCoordinates.y >= 0 && barycentricCoordinates.z >=0)
    {
        // Update barycentric coordinates.
        barycentricCoordinates = restBasisT * (pos - restBasisAncor);
        // Evaluate the implicit functions.
        dist = abs(getValue(pos));
        if (dist < epsilon)
        {
            out_depth = original_depth(gl_FragCoord.z) + travelled;
            out_barycentricWeights = vec4(barycentricCoordinates.xyz, 1 - barycentricCoordinates.x - barycentricCoordinates.y - barycentricCoordinates.z);
            out_tetrahedronID = gl_PrimitiveID;
            out_normal = getGradient(pos); 
            out_pos = dofsSummit[0] * out_barycentricWeights.x + dofsSummit[1] * out_barycentricWeights.y + dofsSummit[2] * out_barycentricWeights.z + dofsSummit[3] * out_barycentricWeights.w;
            return;
        }
		// Step forwards.
		travelled += dist;
		pos += dist * rayDirection;
    }
    discard;
}