/**
 * @file lscm.cpp (COMPLETE IMPLEMENTATION)
 * @brief LSCM (Least Squares Conformal Maps) parameterization
 *
 * Algorithm:
 * 1. Build local vertex mapping (global → local indices)
 * 2. Assemble LSCM sparse matrix
 * 3. Set boundary conditions (pin 2 vertices)
 * 4. Solve sparse linear system
 * 5. Normalize UVs to [0,1]²
 */

#include "lscm.h"
#include "math_utils.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <map>
#include <vector>
#include <set>

// Eigen library for sparse matrices (using third_party/eigen)
#include "../third_party/eigen/Eigen/Sparse"
#include "../third_party/eigen/Eigen/SparseLU"

/**
 * @brief Find boundary vertices in an island
 * 
 * Boundary edges appear only once in the island's face set.
 */
int find_boundary_vertices(const Mesh* mesh,
                          const int* face_indices,
                          int num_faces,
                          int** boundary_out) {
    // Count edge occurrences
    std::map<std::pair<int, int>, int> edge_count;
    
    for (int i = 0; i < num_faces; i++) {
        int face_idx = face_indices[i];
        int v0 = mesh->triangles[3 * face_idx + 0];
        int v1 = mesh->triangles[3 * face_idx + 1];
        int v2 = mesh->triangles[3 * face_idx + 2];
        
        // Add three edges (ensure consistent ordering)
        auto edge1 = std::make_pair(std::min(v0, v1), std::max(v0, v1));
        auto edge2 = std::make_pair(std::min(v1, v2), std::max(v1, v2));
        auto edge3 = std::make_pair(std::min(v2, v0), std::max(v2, v0));
        
        edge_count[edge1]++;
        edge_count[edge2]++;
        edge_count[edge3]++;
    }
    
    // Collect boundary vertices (from edges that appear only once)
    std::set<int> boundary_verts;
    
    for (const auto& pair : edge_count) {
        if (pair.second == 1) {
            boundary_verts.insert(pair.first.first);
            boundary_verts.insert(pair.first.second);
        }
    }
    
    // Convert to array
    int num_boundary = boundary_verts.size();
    *boundary_out = (int*)malloc(num_boundary * sizeof(int));

    int idx = 0;
    for (int v : boundary_verts) {
        (*boundary_out)[idx++] = v;
    }

    // printf("  Found %d boundary vertices\n", num_boundary);
    return num_boundary;
}

void normalize_uvs_to_unit_square(float* uvs, int num_verts) {
    if (!uvs || num_verts == 0) return;

    // Find bounding box
    float min_u = FLT_MAX, max_u = -FLT_MAX;
    float min_v = FLT_MAX, max_v = -FLT_MAX;

    for (int i = 0; i < num_verts; i++) {
        float u = uvs[i * 2];
        float v = uvs[i * 2 + 1];

        min_u = min_float(min_u, u);
        max_u = max_float(max_u, u);
        min_v = min_float(min_v, v);
        max_v = max_float(max_v, v);
    }

    float u_range = max_u - min_u;
    float v_range = max_v - min_v;

    if (u_range < 1e-6f) u_range = 1.0f;
    if (v_range < 1e-6f) v_range = 1.0f;

    // Normalize to [0, 1]
    for (int i = 0; i < num_verts; i++) {
        uvs[i * 2] = (uvs[i * 2] - min_u) / u_range;
        uvs[i * 2 + 1] = (uvs[i * 2 + 1] - min_v) / v_range;
    }
}

/**
 * @brief Add LSCM contribution for a single triangle
 * 
 * @param triplets Matrix entries to add
 * @param v0, v1, v2 Local vertex indices
 * @param p0, p1, p2 3D positions
 */
static void add_triangle_contribution(
    std::vector<Eigen::Triplet<double>>& triplets,
    int v0, int v1, int v2,
    Vec3 p0, Vec3 p1, Vec3 p2)
{
    // Project triangle to its plane
    Vec3 e1 = vec3_sub(p1, p0);
    Vec3 e2 = vec3_sub(p2, p0);
    
    Vec3 normal = vec3_normalize(vec3_cross(e1, e2));
    Vec3 u_axis = vec3_normalize(e1);
    Vec3 v_axis = vec3_cross(normal, u_axis);
    
    // Local 2D coordinates
    double q0_x = 0.0;
    double q0_y = 0.0;
    
    double q1_x = vec3_dot(e1, u_axis);
    double q1_y = vec3_dot(e1, v_axis);
    
    double q2_x = vec3_dot(e2, u_axis);
    double q2_y = vec3_dot(e2, v_axis);
    
    // Triangle area (weight)
    double area = 0.5 * fabs(q1_x * q2_y - q1_y * q2_x);
    
    if (area < 1e-10) return; // Degenerate triangle
    
    // Add LSCM energy terms for each edge
    
    // Edge v0 → v1
    {
        double dx = q1_x - q0_x;
        double dy = q1_y - q0_y;
        
        triplets.push_back(Eigen::Triplet<double>(2*v0,   2*v1,    area * dx));
        triplets.push_back(Eigen::Triplet<double>(2*v0,   2*v1+1,  area * dy));
        triplets.push_back(Eigen::Triplet<double>(2*v0+1, 2*v1,    area * dy));
        triplets.push_back(Eigen::Triplet<double>(2*v0+1, 2*v1+1,  area * (-dx)));
        
        triplets.push_back(Eigen::Triplet<double>(2*v0,   2*v0,    -area * dx));
        triplets.push_back(Eigen::Triplet<double>(2*v0,   2*v0+1,  -area * dy));
        triplets.push_back(Eigen::Triplet<double>(2*v0+1, 2*v0,    -area * dy));
        triplets.push_back(Eigen::Triplet<double>(2*v0+1, 2*v0+1,  -area * (-dx)));
    }
    
    // Edge v1 → v2
    {
        double dx = q2_x - q1_x;
        double dy = q2_y - q1_y;
        
        triplets.push_back(Eigen::Triplet<double>(2*v1,   2*v2,    area * dx));
        triplets.push_back(Eigen::Triplet<double>(2*v1,   2*v2+1,  area * dy));
        triplets.push_back(Eigen::Triplet<double>(2*v1+1, 2*v2,    area * dy));
        triplets.push_back(Eigen::Triplet<double>(2*v1+1, 2*v2+1,  area * (-dx)));
        
        triplets.push_back(Eigen::Triplet<double>(2*v1,   2*v1,    -area * dx));
        triplets.push_back(Eigen::Triplet<double>(2*v1,   2*v1+1,  -area * dy));
        triplets.push_back(Eigen::Triplet<double>(2*v1+1, 2*v1,    -area * dy));
        triplets.push_back(Eigen::Triplet<double>(2*v1+1, 2*v1+1,  -area * (-dx)));
    }
    
    // Edge v2 → v0
    {
        double dx = q0_x - q2_x;
        double dy = q0_y - q2_y;
        
        triplets.push_back(Eigen::Triplet<double>(2*v2,   2*v0,    area * dx));
        triplets.push_back(Eigen::Triplet<double>(2*v2,   2*v0+1,  area * dy));
        triplets.push_back(Eigen::Triplet<double>(2*v2+1, 2*v0,    area * dy));
        triplets.push_back(Eigen::Triplet<double>(2*v2+1, 2*v0+1,  area * (-dx)));
        
        triplets.push_back(Eigen::Triplet<double>(2*v2,   2*v2,    -area * dx));
        triplets.push_back(Eigen::Triplet<double>(2*v2,   2*v2+1,  -area * dy));
        triplets.push_back(Eigen::Triplet<double>(2*v2+1, 2*v2,    -area * dy));
        triplets.push_back(Eigen::Triplet<double>(2*v2+1, 2*v2+1,  -area * (-dx)));
    }
}

/**
 * @brief Find two vertices far apart for pinning
 * 
 * @param mesh Input mesh
 * @param vertices Vertex indices to search
 * @param num_verts Number of vertices
 * @param pin1_out First pin vertex
 * @param pin2_out Second pin vertex
 */
static void find_pin_vertices(const Mesh* mesh,
                               const int* vertices,
                               int num_verts,
                               int* pin1_out,
                               int* pin2_out) {
    if (num_verts < 2) {
        *pin1_out = vertices[0];
        *pin2_out = vertices[0];
        return;
    }
    
    // Find two vertices with maximum distance
    float max_dist = -1.0f;
    int best_i = 0, best_j = 1;
    
    for (int i = 0; i < num_verts; i++) {
        Vec3 pi = get_vertex_position(mesh, vertices[i]);
        
        for (int j = i + 1; j < num_verts; j++) {
            Vec3 pj = get_vertex_position(mesh, vertices[j]);
            Vec3 diff = vec3_sub(pj, pi);
            float dist = vec3_length(diff);
            
            if (dist > max_dist) {
                max_dist = dist;
                best_i = i;
                best_j = j;
            }
        }
    }
    
    *pin1_out = vertices[best_i];
    *pin2_out = vertices[best_j];
}

float* lscm_parameterize(const Mesh* mesh,
                         const int* face_indices,
                         int num_faces) {
    if (!mesh || !face_indices || num_faces == 0) return NULL;

    // printf("LSCM parameterizing %d faces...\n", num_faces);

    // STEP 1: Build local vertex mapping
    std::map<int, int> global_to_local;
    std::vector<int> local_to_global;
    
    for (int i = 0; i < num_faces; i++) {
        int face_idx = face_indices[i];
        int v0 = mesh->triangles[3 * face_idx + 0];
        int v1 = mesh->triangles[3 * face_idx + 1];
        int v2 = mesh->triangles[3 * face_idx + 2];
        
        if (global_to_local.find(v0) == global_to_local.end()) {
            global_to_local[v0] = local_to_global.size();
            local_to_global.push_back(v0);
        }
        if (global_to_local.find(v1) == global_to_local.end()) {
            global_to_local[v1] = local_to_global.size();
            local_to_global.push_back(v1);
        }
        if (global_to_local.find(v2) == global_to_local.end()) {
            global_to_local[v2] = local_to_global.size();
            local_to_global.push_back(v2);
        }
    }
    
    int n = local_to_global.size();
    // printf("  Island has %d vertices\n", n);

    if (n < 3) {
        // fprintf(stderr, "LSCM: Island too small (%d vertices)\n", n);
        return NULL;
    }

    // STEP 2: Build LSCM sparse matrix
    typedef Eigen::Triplet<double> T;
    std::vector<T> triplets;
    triplets.reserve(num_faces * 72); // 3 edges * 24 entries per edge
    
    for (int i = 0; i < num_faces; i++) {
        int face_idx = face_indices[i];
        int v0_global = mesh->triangles[3 * face_idx + 0];
        int v1_global = mesh->triangles[3 * face_idx + 1];
        int v2_global = mesh->triangles[3 * face_idx + 2];
        
        int v0_local = global_to_local[v0_global];
        int v1_local = global_to_local[v1_global];
        int v2_local = global_to_local[v2_global];
        
        Vec3 p0 = get_vertex_position(mesh, v0_global);
        Vec3 p1 = get_vertex_position(mesh, v1_global);
        Vec3 p2 = get_vertex_position(mesh, v2_global);
        
        add_triangle_contribution(triplets, v0_local, v1_local, v2_local, p0, p1, p2);
    }
    
    // printf("  Built matrix with %zu entries\n", triplets.size());

    // STEP 3: Find vertices to pin
    int* boundary_verts = NULL;
    int num_boundary = find_boundary_vertices(mesh, face_indices, num_faces, &boundary_verts);
    
    int pin1_global, pin2_global;
    
    if (num_boundary >= 2) {
        // Use boundary vertices
        find_pin_vertices(mesh, boundary_verts, num_boundary, &pin1_global, &pin2_global);
    } else {
        // Closed mesh or single boundary - pick any two far apart
        find_pin_vertices(mesh, local_to_global.data(), n, &pin1_global, &pin2_global);
    }
    
    if (boundary_verts) free(boundary_verts);
    
    int pin1_local = global_to_local[pin1_global];
    int pin2_local = global_to_local[pin2_global];
    
    // printf("  Pinning vertices %d (local %d) and %d (local %d)\n", 
    //        pin1_global, pin1_local, pin2_global, pin2_local);

    // STEP 4: Set boundary conditions
    // Clear rows for pinned vertices and set constraints
    
    // Remove all entries in pin1 and pin2 rows
    std::vector<T> filtered_triplets;
    filtered_triplets.reserve(triplets.size());
    
    for (const auto& t : triplets) {
        int row = t.row();
        if (row != 2*pin1_local && row != 2*pin1_local+1 &&
            row != 2*pin2_local && row != 2*pin2_local+1) {
            filtered_triplets.push_back(t);
        }
    }
    
    // Add identity for pinned vertices
    // Pin vertex 1 to (0, 0)
    filtered_triplets.push_back(T(2*pin1_local, 2*pin1_local, 1.0));
    filtered_triplets.push_back(T(2*pin1_local+1, 2*pin1_local+1, 1.0));
    
    // Pin vertex 2 to (1, 0)
    filtered_triplets.push_back(T(2*pin2_local, 2*pin2_local, 1.0));
    filtered_triplets.push_back(T(2*pin2_local+1, 2*pin2_local+1, 1.0));

    // Build matrix and RHS
    Eigen::SparseMatrix<double> A(2*n, 2*n);
    A.setFromTriplets(filtered_triplets.begin(), filtered_triplets.end());

    Eigen::VectorXd b = Eigen::VectorXd::Zero(2*n);
    b[2*pin1_local] = 0.0;     // u = 0
    b[2*pin1_local+1] = 0.0;   // v = 0
    b[2*pin2_local] = 1.0;     // u = 1
    b[2*pin2_local+1] = 0.0;   // v = 0

    // STEP 5: Solve linear system
    // printf("  Solving %dx%d sparse system...\n", 2*n, 2*n);
    
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(A);
    
    if (solver.info() != Eigen::Success) {
        // fprintf(stderr, "LSCM: Matrix decomposition failed\n");
        return NULL;
    }
    
    Eigen::VectorXd x = solver.solve(b);
    
    if (solver.info() != Eigen::Success) {
        // fprintf(stderr, "LSCM: Solve failed\n");
        return NULL;
    }
    
    // printf("  Solve completed\n");

    // STEP 6: Extract UVs
    float* uvs = (float*)malloc(n * 2 * sizeof(float));
    
    for (int i = 0; i < n; i++) {
        uvs[i * 2] = static_cast<float>(x[2*i]);
        uvs[i * 2 + 1] = static_cast<float>(x[2*i + 1]);
    }

    // STEP 7: Normalize to [0,1]²
    normalize_uvs_to_unit_square(uvs, n);

    // printf("  LSCM completed successfully\n");
    return uvs;
}
