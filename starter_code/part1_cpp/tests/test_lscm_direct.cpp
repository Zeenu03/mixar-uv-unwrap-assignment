/**
 * @file test_lscm_direct.cpp
 * @brief Direct test for LSCM parameterization
 * 
 * This test directly calls lscm_parameterize on a simple mesh
 * to verify the implementation is working.
 */

#include "mesh.h"
#include "topology.h"
#include "lscm.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define TEST_DATA_DIR "../../test_data/meshes/"

void test_lscm_direct(const char* mesh_name) {
    printf("\n========================================\n");
    printf("Testing LSCM on: %s\n", mesh_name);
    printf("========================================\n");

    char filename[256];
    snprintf(filename, sizeof(filename), "%s%s", TEST_DATA_DIR, mesh_name);

    // Load mesh
    Mesh* mesh = load_obj(filename);
    if (!mesh) {
        printf("FAIL: Could not load mesh\n");
        return;
    }

    printf("Loaded mesh: %d vertices, %d triangles\n", 
           mesh->num_vertices, mesh->num_triangles);

    // Create face indices array (all faces)
    int* face_indices = (int*)malloc(mesh->num_triangles * sizeof(int));
    for (int i = 0; i < mesh->num_triangles; i++) {
        face_indices[i] = i;
    }

    // Call LSCM
    printf("\nCalling LSCM parameterization...\n");
    float* uvs = lscm_parameterize(mesh, face_indices, mesh->num_triangles);

    if (!uvs) {
        printf("FAIL: LSCM returned NULL\n");
        free(face_indices);
        free_mesh(mesh);
        return;
    }

    printf("\n✓ LSCM succeeded!\n");

    // Analyze results
    printf("\nUV Coordinates Analysis:\n");
    
    float min_u = 1e10f, max_u = -1e10f;
    float min_v = 1e10f, max_v = -1e10f;
    
    for (int i = 0; i < mesh->num_vertices; i++) {
        float u = uvs[i * 2];
        float v = uvs[i * 2 + 1];
        
        if (u < min_u) min_u = u;
        if (u > max_u) max_u = u;
        if (v < min_v) min_v = v;
        if (v > max_v) max_v = v;
    }
    
    printf("  U range: [%.3f, %.3f]\n", min_u, max_u);
    printf("  V range: [%.3f, %.3f]\n", min_v, max_v);
    
    // Check if normalized to [0,1]
    bool properly_normalized = (min_u >= -0.01f && max_u <= 1.01f && 
                                min_v >= -0.01f && max_v <= 1.01f);
    
    if (properly_normalized) {
        printf("  ✓ UVs properly normalized to [0,1]²\n");
    } else {
        printf("  ✗ WARNING: UVs not in [0,1]² range!\n");
    }
    
    // Sample some UV coordinates
    printf("\nSample UV coordinates (first 5 vertices):\n");
    for (int i = 0; i < mesh->num_vertices && i < 5; i++) {
        printf("  v%d: (%.3f, %.3f)\n", i, uvs[i*2], uvs[i*2+1]);
    }
    
    // Check for NaN or infinite values
    int bad_count = 0;
    for (int i = 0; i < mesh->num_vertices; i++) {
        float u = uvs[i * 2];
        float v = uvs[i * 2 + 1];
        
        if (isnan(u) || isnan(v) || isinf(u) || isinf(v)) {
            bad_count++;
        }
    }
    
    if (bad_count > 0) {
        printf("  ✗ WARNING: Found %d vertices with NaN/Inf values!\n", bad_count);
    } else {
        printf("  ✓ All UV values are valid (no NaN/Inf)\n");
    }
    
    // Check for zero-area (all same coordinates)
    bool all_same = true;
    float first_u = uvs[0];
    float first_v = uvs[1];
    for (int i = 1; i < mesh->num_vertices; i++) {
        if (fabs(uvs[i*2] - first_u) > 1e-6f || fabs(uvs[i*2+1] - first_v) > 1e-6f) {
            all_same = false;
            break;
        }
    }
    
    if (all_same) {
        printf("  ✗ WARNING: All UVs are the same (degenerate parameterization)!\n");
    } else {
        printf("  ✓ UVs have variation (non-degenerate)\n");
    }

    printf("\n========================================\n");
    printf("LSCM Test Result: %s\n", 
           (properly_normalized && bad_count == 0 && !all_same) ? "PASS ✓" : "PARTIAL/FAIL ✗");
    printf("========================================\n\n");

    // Cleanup
    free(uvs);
    free(face_indices);
    free_mesh(mesh);
}

int main() {
    printf("\n");
    printf("╔════════════════════════════════════════╗\n");
    printf("║   LSCM Direct Parameterization Test   ║\n");
    printf("╔════════════════════════════════════════╗\n");

    // Test on different meshes
    test_lscm_direct("01_cube.obj");
    test_lscm_direct("02_cylinder.obj");
    test_lscm_direct("03_sphere.obj");

    return 0;
}
