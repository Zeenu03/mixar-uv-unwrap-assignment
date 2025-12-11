/**
 * @file benchmark_lscm_performance.cpp
 * @brief Performance benchmark for LSCM parameterization
 * 
 * Performance Targets:
 * - 1k vertices with SparseLU: < 0.5s
 * - 10k vertices with SparseLU: < 5s
 */

#include "mesh.h"
#include "lscm.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

/**
 * Generate a planar grid mesh for benchmarking
 */
Mesh* generate_grid_mesh(int grid_size) {
    int num_verts = grid_size * grid_size;
    int num_tris = 2 * (grid_size - 1) * (grid_size - 1);
    
    Mesh* mesh = (Mesh*)malloc(sizeof(Mesh));
    mesh->num_vertices = num_verts;
    mesh->num_triangles = num_tris;
    mesh->vertices = (float*)malloc(3 * num_verts * sizeof(float));
    mesh->triangles = (int*)malloc(3 * num_tris * sizeof(int));
    mesh->uvs = NULL;
    
    // Generate grid vertices
    float spacing = 1.0f / (grid_size - 1);
    for (int i = 0; i < grid_size; i++) {
        for (int j = 0; j < grid_size; j++) {
            int idx = i * grid_size + j;
            mesh->vertices[3*idx + 0] = j * spacing;
            mesh->vertices[3*idx + 1] = i * spacing;
            mesh->vertices[3*idx + 2] = 0.0f;
        }
    }
    
    // Generate triangles
    int tri_idx = 0;
    for (int i = 0; i < grid_size - 1; i++) {
        for (int j = 0; j < grid_size - 1; j++) {
            int v0 = i * grid_size + j;
            int v1 = i * grid_size + (j + 1);
            int v2 = (i + 1) * grid_size + j;
            int v3 = (i + 1) * grid_size + (j + 1);
            
            mesh->triangles[3*tri_idx + 0] = v0;
            mesh->triangles[3*tri_idx + 1] = v1;
            mesh->triangles[3*tri_idx + 2] = v2;
            tri_idx++;
            
            mesh->triangles[3*tri_idx + 0] = v1;
            mesh->triangles[3*tri_idx + 1] = v3;
            mesh->triangles[3*tri_idx + 2] = v2;
            tri_idx++;
        }
    }
    
    return mesh;
}

double get_time_in_seconds() {
    return (double)clock() / CLOCKS_PER_SEC;
}

int main() {
    printf("\n");
    printf("========================================================\n");
    printf("         LSCM Performance Benchmark Test\n");
    printf("  Targets: 1k verts < 0.5s, 10k verts < 5s\n");
    printf("========================================================\n");
    
    int passed = 0;
    int failed = 0;
    
    // Test 1: ~1k vertices (32x32 grid = 1024 vertices)
    printf("\n[Test 1/2] 1k vertices (32x32 grid = 1024 vertices)\n");
    printf("Target: < 0.5 seconds\n");
    printf("--------------------------------------------------------\n");
    
    Mesh* mesh_1k = generate_grid_mesh(32);
    int* faces_1k = (int*)malloc(mesh_1k->num_triangles * sizeof(int));
    for (int i = 0; i < mesh_1k->num_triangles; i++) faces_1k[i] = i;
    
    double start_1k = get_time_in_seconds();
    float* uvs_1k = lscm_parameterize(mesh_1k, faces_1k, mesh_1k->num_triangles);
    double time_1k = get_time_in_seconds() - start_1k;
    
    if (uvs_1k) {
        printf("Time elapsed: %.4f seconds\n", time_1k);
        if (time_1k <= 0.5) {
            printf("PASS - Within target (%.2fx faster than target)\n", 0.5 / time_1k);
            passed++;
        } else {
            printf("FAIL - Exceeded target by %.3f seconds\n", time_1k - 0.5);
            failed++;
        }
        free(uvs_1k);
    } else {
        printf("FAIL - LSCM failed\n");
        failed++;
    }
    
    free(faces_1k);
    free(mesh_1k->vertices);
    free(mesh_1k->triangles);
    free(mesh_1k);
    
    // Test 2: ~10k vertices (100x100 grid = 10000 vertices)
    printf("\n[Test 2/2] 10k vertices (100x100 grid = 10000 vertices)\n");
    printf("Target: < 5 seconds\n");
    printf("--------------------------------------------------------\n");
    
    Mesh* mesh_10k = generate_grid_mesh(100);
    int* faces_10k = (int*)malloc(mesh_10k->num_triangles * sizeof(int));
    for (int i = 0; i < mesh_10k->num_triangles; i++) faces_10k[i] = i;
    
    double start_10k = get_time_in_seconds();
    float* uvs_10k = lscm_parameterize(mesh_10k, faces_10k, mesh_10k->num_triangles);
    double time_10k = get_time_in_seconds() - start_10k;
    
    if (uvs_10k) {
        printf("Time elapsed: %.4f seconds\n", time_10k);
        if (time_10k <= 5.0) {
            printf("PASS - Within target (%.2fx faster than target)\n", 5.0 / time_10k);
            passed++;
        } else {
            printf("FAIL - Exceeded target by %.3f seconds\n", time_10k - 5.0);
            failed++;
        }
        free(uvs_10k);
    } else {
        printf("FAIL - LSCM failed\n");
        failed++;
    }
    
    free(faces_10k);
    free(mesh_10k->vertices);
    free(mesh_10k->triangles);
    free(mesh_10k);
    
    // Summary
    printf("\n========================================================\n");
    printf("              Performance Summary\n");
    printf("========================================================\n");
    printf("  Tests Passed: %d/2\n", passed);
    printf("  Tests Failed: %d/2\n", failed);
    printf("========================================================\n");
    
    if (failed == 0) {
        printf("\nAll performance targets MET!\n\n");
    } else {
        printf("\nSome performance targets NOT MET\n\n");
    }
    
    return (failed == 0) ? 0 : 1;
}
