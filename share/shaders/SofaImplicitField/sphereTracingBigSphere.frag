#version 460
// Fragment shadder.

// Replaces "sphereTracing_v1.frag" in the pipeline:
// sphereTracing_v1.vert -> sphereTracing_v1.geom -> sphereTracing_v1.frag

// Hard coded implicit function of a sphere with center (3,0,0) and radius 4.

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

uniform vec3 viewingDirection;
uniform float epsilon;
uniform float h;
uniform float zNear;
uniform float zFar;

float unormalized_depth(float d)
{
    return d * abs(zFar - zNear);
}

float normalized_depth(float d)
{
    return (d - zNear) / (zFar - zNear);
}

float getValue(in vec3 pos)
{
    vec3 center = {3, 0, 0};
    return length(pos-center) - 4;
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
    // Begin sphere tracing.
    vec3 pos = rayOrigin;
    vec3 barycentricCoordinates = {0,0,0};
    float dist;

    // Ray march for as long as the current position is inside the target domain.
    while (barycentricCoordinates.x + barycentricCoordinates.y + barycentricCoordinates.z <= 1 && barycentricCoordinates.x >= 0 && barycentricCoordinates.y >= 0 && barycentricCoordinates.z >=0)
    {
        // update barycentric coordinates.
        barycentricCoordinates = restBasisT * (pos - restBasisAncor);
        // Evaluate the implicit functions.
        dist = abs(getValue(pos));
        if (dist < epsilon)
        {
            out_barycentricWeights = vec4(barycentricCoordinates.xyz, 1 - barycentricCoordinates.x - barycentricCoordinates.y - barycentricCoordinates.z);
            out_tetrahedronID = gl_PrimitiveID;
            out_normal = getGradient(pos); 
            out_pos = dofsSummit[0] * out_barycentricWeights.x + dofsSummit[1] * out_barycentricWeights.y + dofsSummit[2] * out_barycentricWeights.z + dofsSummit[3] * out_barycentricWeights.w;
            float d = normalized_depth(dot(out_pos, abs(viewingDirection)));
            out_depth = unormalized_depth(d);
            gl_FragDepth = d;

            return;
        }
		// Step forwards.
		pos += dist * rayDirection;
    }
    discard;
}