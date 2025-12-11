/**
 * @file packing.cpp (COMPLETE IMPLEMENTATION)
 * @brief UV island packing into [0,1]² texture space
 *
 * Algorithm: Shelf packing
 * 1. Compute bounding box for each island
 * 2. Sort islands by height (descending)
 * 3. Pack using shelf algorithm
 * 4. Scale to fit [0,1]²
 * 
 * Expected coverage: > 60%
 */

#include "unwrap.h"
#include "math_utils.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <vector>
#include <algorithm>

/**
 * @brief Island bounding box info
 */
struct Island {
    int id;
    float min_u, max_u, min_v, max_v;
    float width, height;
    float target_x, target_y;  // Packed position
    std::vector<int> vertex_indices;
};

/**
 * @brief Shelf structure for packing
 */
struct Shelf {
    float x;           // Current x position in shelf
    float y;           // Y position of shelf
    float height;      // Height of this shelf
    float max_width;   // Maximum width (usually 1.0 before scaling)
    
    Shelf(float y_pos, float h) : x(0.0f), y(y_pos), height(h), max_width(10000.0f) {}
};

void pack_uv_islands(Mesh* mesh,
                     const UnwrapResult* result,
                     float margin) {
    if (!mesh || !result || !mesh->uvs) return;

    if (result->num_islands <= 1) {
        // Single island, already normalized to [0,1]
        printf("Single island, no packing needed\n");
        return;
    }

    // printf("Packing %d islands with margin %.3f...\n", result->num_islands, margin);

    std::vector<Island> islands(result->num_islands);

    // STEP 1: Compute bounding box for each island
    // printf("  Computing island bounding boxes...\n");
    
    // Initialize islands
    for (int i = 0; i < result->num_islands; i++) {
        islands[i].id = i;
        islands[i].min_u = FLT_MAX;
        islands[i].max_u = -FLT_MAX;
        islands[i].min_v = FLT_MAX;
        islands[i].max_v = -FLT_MAX;
    }
    
    // Find vertices in each island and compute bounds
    for (int face_idx = 0; face_idx < mesh->num_triangles; face_idx++) {
        int island_id = result->face_island_ids[face_idx];
        
        // Get face vertices
        int v0 = mesh->triangles[3 * face_idx + 0];
        int v1 = mesh->triangles[3 * face_idx + 1];
        int v2 = mesh->triangles[3 * face_idx + 2];
        
        int verts[3] = {v0, v1, v2};
        
        for (int i = 0; i < 3; i++) {
            int v = verts[i];
            float u = mesh->uvs[2 * v + 0];
            float v_coord = mesh->uvs[2 * v + 1];
            
            // Update bounds
            islands[island_id].min_u = min_float(islands[island_id].min_u, u);
            islands[island_id].max_u = max_float(islands[island_id].max_u, u);
            islands[island_id].min_v = min_float(islands[island_id].min_v, v_coord);
            islands[island_id].max_v = max_float(islands[island_id].max_v, v_coord);
            
            // Track vertex belongs to this island
            bool already_added = false;
            for (int vi : islands[island_id].vertex_indices) {
                if (vi == v) {
                    already_added = true;
                    break;
                }
            }
            if (!already_added) {
                islands[island_id].vertex_indices.push_back(v);
            }
        }
    }
    
    // Compute dimensions
    for (int i = 0; i < result->num_islands; i++) {
        islands[i].width = islands[i].max_u - islands[i].min_u;
        islands[i].height = islands[i].max_v - islands[i].min_v;
        printf("  Island %d: %.3f x %.3f (%d vertices)\n", 
               i, islands[i].width, islands[i].height, 
               (int)islands[i].vertex_indices.size());
    }

    // STEP 2: Sort islands by height (descending - tallest first)
    // printf("  Sorting islands by height...\n");
    std::sort(islands.begin(), islands.end(), 
              [](const Island& a, const Island& b) {
                  return a.height > b.height;
              });

    // STEP 3: Shelf packing
    // printf("  Packing islands using shelf algorithm...\n");
    
    std::vector<Shelf> shelves;
    shelves.push_back(Shelf(0.0f, islands[0].height));
    
    for (int i = 0; i < islands.size(); i++) {
        Island& island = islands[i];
        Shelf& current_shelf = shelves.back();
        
        // Check if island fits in current shelf
        if (current_shelf.x + island.width + margin <= current_shelf.max_width) {
            // Fits in current shelf
            island.target_x = current_shelf.x;
            island.target_y = current_shelf.y;
            current_shelf.x += island.width + margin;
        } else {
            // Need new shelf
            float new_y = current_shelf.y + current_shelf.height + margin;
            shelves.push_back(Shelf(new_y, island.height));
            
            island.target_x = 0.0f;
            island.target_y = new_y;
            shelves.back().x = island.width + margin;
        }
        
        // printf("    Island %d placed at (%.3f, %.3f)\n", 
            //    island.id, island.target_x, island.target_y);
    }
    
    // printf("  Used %zu shelves\n", shelves.size());

    // STEP 4: Move islands to packed positions
    // printf("  Moving islands to packed positions...\n");
    
    for (const Island& island : islands) {
        float offset_x = island.target_x - island.min_u;
        float offset_y = island.target_y - island.min_v;
        
        // Apply offset to all vertices in this island
        for (int v : island.vertex_indices) {
            mesh->uvs[2 * v + 0] += offset_x;
            mesh->uvs[2 * v + 1] += offset_y;
        }
    }

    // STEP 5: Scale to fit [0,1]²
    // printf("  Scaling to [0,1]²...\n");
    
    // Find bounding box of entire packed result
    float global_min_u = FLT_MAX, global_max_u = -FLT_MAX;
    float global_min_v = FLT_MAX, global_max_v = -FLT_MAX;
    
    for (int v = 0; v < mesh->num_vertices; v++) {
        float u = mesh->uvs[2 * v + 0];
        float v_coord = mesh->uvs[2 * v + 1];
        
        global_min_u = min_float(global_min_u, u);
        global_max_u = max_float(global_max_u, u);
        global_min_v = min_float(global_min_v, v_coord);
        global_max_v = max_float(global_max_v, v_coord);
    }
    
    float packed_width = global_max_u - global_min_u;
    float packed_height = global_max_v - global_min_v;
    float scale = 1.0f / max_float(packed_width, packed_height);
    
    // printf("  Packed dimensions: %.3f x %.3f\n", packed_width, packed_height);
    // printf("  Scale factor: %.3f\n", scale);
    
    // Apply scaling
    for (int v = 0; v < mesh->num_vertices; v++) {
        mesh->uvs[2 * v + 0] = (mesh->uvs[2 * v + 0] - global_min_u) * scale;
        mesh->uvs[2 * v + 1] = (mesh->uvs[2 * v + 1] - global_min_v) * scale;
    }

    // printf("  Packing completed successfully\n");
}

void compute_quality_metrics(const Mesh* mesh, UnwrapResult* result) {
    if (!mesh || !result || !mesh->uvs) return;

    // printf("Computing quality metrics...\n");

    // COVERAGE: Approximate using island bounding boxes
    float total_area = 0.0f;
    
    for (int island_id = 0; island_id < result->num_islands; island_id++) {
        float min_u = FLT_MAX, max_u = -FLT_MAX;
        float min_v = FLT_MAX, max_v = -FLT_MAX;
        
        // Find bounds of this island
        for (int face_idx = 0; face_idx < mesh->num_triangles; face_idx++) {
            if (result->face_island_ids[face_idx] == island_id) {
                int v0 = mesh->triangles[3 * face_idx + 0];
                int v1 = mesh->triangles[3 * face_idx + 1];
                int v2 = mesh->triangles[3 * face_idx + 2];
                
                for (int v : {v0, v1, v2}) {
                    float u = mesh->uvs[2 * v + 0];
                    float v_coord = mesh->uvs[2 * v + 1];
                    
                    min_u = min_float(min_u, u);
                    max_u = max_float(max_u, u);
                    min_v = min_float(min_v, v_coord);
                    max_v = max_float(max_v, v_coord);
                }
            }
        }
        
        float width = max_u - min_u;
        float height = max_v - min_v;
        total_area += width * height;
    }
    
    result->coverage = total_area; // Approximate (actual would need rasterization)
    
    // STRETCH: Simple estimation (would need SVD for accurate computation)
    float sum_stretch = 0.0f;
    float max_stretch = 0.0f;
    int count = 0;
    
    for (int face_idx = 0; face_idx < mesh->num_triangles; face_idx++) {
        int v0 = mesh->triangles[3 * face_idx + 0];
        int v1 = mesh->triangles[3 * face_idx + 1];
        int v2 = mesh->triangles[3 * face_idx + 2];
        
        // Get 3D positions
        Vec3 p0 = get_vertex_position(mesh, v0);
        Vec3 p1 = get_vertex_position(mesh, v1);
        Vec3 p2 = get_vertex_position(mesh, v2);
        
        // Get UV coordinates
        float u0 = mesh->uvs[2 * v0 + 0];
        float v0_coord = mesh->uvs[2 * v0 + 1];
        float u1 = mesh->uvs[2 * v1 + 0];
        float v1_coord = mesh->uvs[2 * v1 + 1];
        float u2 = mesh->uvs[2 * v2 + 0];
        float v2_coord = mesh->uvs[2 * v2 + 1];
        
        // Compute edge lengths in 3D
        Vec3 e1_3d = vec3_sub(p1, p0);
        Vec3 e2_3d = vec3_sub(p2, p0);
        float len1_3d = vec3_length(e1_3d);
        float len2_3d = vec3_length(e2_3d);
        
        // Compute edge lengths in UV
        float du1 = u1 - u0;
        float dv1 = v1_coord - v0_coord;
        float du2 = u2 - u0;
        float dv2 = v2_coord - v0_coord;
        
        float len1_uv = sqrtf(du1*du1 + dv1*dv1);
        float len2_uv = sqrtf(du2*du2 + dv2*dv2);
        
        // Simple stretch estimate (ratio of UV to 3D lengths)
        if (len1_3d > 1e-6f && len2_3d > 1e-6f) {
            float stretch1 = len1_uv / len1_3d;
            float stretch2 = len2_uv / len2_3d;
            float tri_stretch = max_float(stretch1, stretch2) / min_float(stretch1, stretch2);
            
            sum_stretch += tri_stretch;
            max_stretch = max_float(max_stretch, tri_stretch);
            count++;
        }
    }
    
    result->avg_stretch = (count > 0) ? sum_stretch / count : 1.0f;
    result->max_stretch = max_stretch;

    // printf("  Avg stretch: %.2f\n", result->avg_stretch);
    // printf("  Max stretch: %.2f\n", result->max_stretch);
    // printf("  Coverage: %.1f%%\n", result->coverage * 100);
}
