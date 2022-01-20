#version 460
// GeometricShadder shadder.

// Belongs to the OpenGL pipeline:
// sphereTracing_v1.vert -> sphereTracing_v1.geom -> sphereTracing_v1.frag

layout (lines_adjacency) in; // This feels like cheating. but allows me to read 4 vertices at the time.
layout (triangle_strip, max_vertices = 6) out;

// in int gl_PrimitiveIDIn;
// in gl_PerVertex {
//  vec4 gl_Position;
//  ...
//  } gl_in[];
in int vertexID[];

// out int gl_PrimitiveID;
out vec3 rayOrigin;
flat out vec3 rayDirection;
flat out mat3 restBasisT;
flat out vec3 restBasisAncor;
flat out vec3[4] dofsSummit;

uniform vec3 viewingDirection;
uniform float smallStep;
uniform samplerBuffer restBuffer; // Vertex rest-position array as a buffer texture.
uniform samplerBuffer dofBuffer; // Vertex position array as a buffer texture.


void main()
{
    // Read ~un-projected~ current and rest-positions of the vertices.
    // ... deformed tetrahedron:
    vec3 posA = texelFetch(dofBuffer, vertexID[0]).xyz;
    vec3 posB = texelFetch(dofBuffer, vertexID[1]).xyz;
    vec3 posC = texelFetch(dofBuffer, vertexID[2]).xyz;
    vec3 posD = texelFetch(dofBuffer, vertexID[3]).xyz;
    // ... undeformed tetrahedron:
    vec3 restA = texelFetch(restBuffer, vertexID[0]).xyz;
    vec3 restB = texelFetch(restBuffer, vertexID[1]).xyz;
    vec3 restC = texelFetch(restBuffer, vertexID[2]).xyz;
    vec3 restD = texelFetch(restBuffer, vertexID[3]).xyz;

    // Compute the basis (mat3 T) of the deformed and undeformed tetrahedron.
    // ... deformed tetrahedron:
    mat3 currentBasisT;
    currentBasisT[0] = posA-posD;
    currentBasisT[1] = posB-posD;
    currentBasisT[2] = posC-posD;
    currentBasisT = inverse(currentBasisT);
    // ... undeformed tetrahedron:
    restBasisT[0] = restA-restD;
    restBasisT[1] = restB-restD;
    restBasisT[2] = restC-restD;
    restBasisT = inverse(restBasisT);

    // Pass through the tetrahedron summits.
    // ... deformed tetrahedron:
    dofsSummit[0] = posA;
    dofsSummit[1] = posB;
    dofsSummit[2] = posC;
    dofsSummit[3] = posD;
    // ... undeformed tetrahedron:
    restBasisAncor = restD;

    // Compute the ray direction in the undeformed tetrahedron.
    vec3 center = 0.25 * posA + 0.25 * posB + 0.25 * posC + 0.25 * posD;
    vec3 center_ = 0.25 * restA + 0.25 * restB + 0.25 * restC + 0.25 * restD;

    vec3 next = center + viewingDirection * smallStep;
    vec3 nextCoefs = currentBasisT * (next - posD);
    vec3 next_ = nextCoefs[0] * restA + nextCoefs[1] * restB + nextCoefs[2] * restC + (1-nextCoefs[0]-nextCoefs[1]-nextCoefs[2]) * restD;

    rayDirection = normalize(next_ - center_);

    // Pass through the ID of the tetrahedron.
    // (( gl_PrimitiveID = gl_PrimitiveIDIn; )) for each vertex before emitting it.
    
    // Draw the tetrahedron's faces using the projected vertex coordinates.
    // Triangle strip: <A,B,C,D,A,B>
    
    // Compute the ray origin through interpolation of undeformed tetrahedron summits.

    // <A,B,C>
    gl_PrimitiveID = gl_PrimitiveIDIn;
    gl_Position = gl_in[0].gl_Position;
    rayOrigin = restA;
    EmitVertex();

    gl_PrimitiveID = gl_PrimitiveIDIn;
    gl_Position = gl_in[1].gl_Position;
    rayOrigin = restB;
    EmitVertex();
    
    gl_PrimitiveID = gl_PrimitiveIDIn;
    gl_Position = gl_in[2].gl_Position;
    rayOrigin = restC;
    EmitVertex();
    
    // <(B,C),D>
    gl_PrimitiveID = gl_PrimitiveIDIn;
    gl_Position = gl_in[3].gl_Position;
    rayOrigin = restD;
    EmitVertex();

    // <(C,D),A>
    gl_PrimitiveID = gl_PrimitiveIDIn;
    gl_Position = gl_in[0].gl_Position;
    rayOrigin = restA;
    EmitVertex();

    // <(D,A),B>
    gl_PrimitiveID = gl_PrimitiveIDIn;
    gl_Position = gl_in[1].gl_Position;
    rayOrigin = restB;
    EmitVertex();

    EndPrimitive();

    return;
}