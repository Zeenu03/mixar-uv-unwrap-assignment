/**
 * @file seam_detection.cpp
 * @brief Seam detection using spanning tree + angular defect
 *
 * SKELETON - YOU IMPLEMENT THIS
 *
 * Algorithm:
 * 1. Build dual graph (faces as nodes, shared edges as edges)
 * 2. Compute spanning tree via BFS
 * 3. Mark non-tree edges as seam candidates
 * 4. Refine using angular defect
 *
 * See reference/algorithms.md for detailed description
 */

#include "unwrap.h"
#include "math_utils.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <set>
#include <queue>
#include <unordered_map>
#include <algorithm>
using namespace std;
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
/**
 * @brief Compute angular defect at a vertex
 *
 * Angular defect = 2π - sum of angles at vertex
 *
 * - Flat surface: defect ≈ 0
 * - Corner (like cube): defect > 0
 * - Saddle: defect < 0
 *
 * @param mesh Input mesh
 * @param vertex_idx Vertex index
 * @return Angular defect in radians
 */
static float compute_angular_defect(const Mesh* mesh, int vertex_idx) {
    // TODO: Implement
    //
    // Steps:
    // 1. Find all triangles containing this vertex
    // 2. For each triangle, compute angle at this vertex
    //    (use compute_vertex_angle_in_triangle from math_utils.h)
    // 3. Sum all angles
    // 4. Return 2*PI - sum
    //
    // Hint: Angular defect indicates curvature
    //       High defect → sharp feature → good seam location

    float angle_sum = 0.0f;

    // YOUR CODE HERE
    for(int face_idx = 0; face_idx < mesh->num_triangles; face_idx++) {
        int v0 = mesh->triangles[3 * face_idx + 0];
        int v1 = mesh->triangles[3 * face_idx + 1];
        int v2 = mesh->triangles[3 * face_idx + 2];
        
        if(v0 == vertex_idx || v1 == vertex_idx || v2 == vertex_idx) {
            float angle = compute_vertex_angle_in_triangle(mesh, face_idx, vertex_idx);
            angle_sum += angle;
        }
    }

    return static_cast<float>(2.0f * M_PI - angle_sum);
}

/**
 * @brief Get all edges incident to a vertex
 */
static std::vector<int> get_vertex_edges(const TopologyInfo* topo, int vertex_idx) {
    std::vector<int> edges;

    // TODO: Iterate through all edges, add those touching vertex_idx

    // YOUR CODE HERE
    for(int idx = 0; idx < topo->num_edges; idx++) {
        int v0 = topo->edges[2 * idx + 0];
        int v1 = topo->edges[2 * idx + 1];

        if(v0 == vertex_idx || v1 == vertex_idx) {
            edges.push_back(idx);
        }
    }

    return edges;
}

long long make_edge_key(int f0, int f1) {
    if(f0 > f1) {
        swap(f0, f1);
    }

    return (static_cast<long long>(f0) << 32) | static_cast<long long>(f1);
}


struct SeamCandidate {
    int edge_idx;
    float score;
};


int find_parent(int v, vector<int>& parent) {
    if(parent[v] == v) {
        return v;
    }
    return parent[v] = find_parent(parent[v], parent);
}

void unite(int a, int b, vector<int>& parent) {
    int pa = find_parent(a, parent);
    int pb = find_parent(b, parent);
    if(pa != pb) {
        parent[pa] = pb;
    }
}

int* detect_seams(const Mesh* mesh,
                  const TopologyInfo* topo,
                  float angle_threshold,
                  int* num_seams_out) {
    if (!mesh || !topo || !num_seams_out) return NULL;
    float angle_threshold_redians =  angle_threshold * static_cast<float>(M_PI)/180.0f;              
    

    // TODO: Implement seam detection
    //
    // ALGORITHM:
    //
    // STEP 1: Build dual graph
    //   - Nodes = faces
    //   - Edges = shared edges between faces
    //   - Use std::vector<std::vector<int>> for adjacency list
    //
    // STEP 2: Spanning tree via BFS
    //   - Start from face 0
    //   - Mark edges in spanning tree
    //   - Use std::set<int> to track tree edges
    //
    // STEP 3: Initial seam candidates = non-tree edges
    //   - All edges NOT in spanning tree
    //
    // STEP 4: Angular defect refinement
    //   - For each vertex with high angular defect (> 0.5 radians)
    //   - Add incident edges to seam candidates
    //
    // STEP 5: Convert seam candidates to array
    //
    // Expected seam counts:
    //   Cube: 7-9 seams
    //   Sphere: 1-3 seams
    //   Cylinder: 1-2 seams

    std::set<int> seam_candidates;
    int boundry_edges_count = 0;
    // STEP 1: Build dual graph
    int F = mesh->num_triangles;
    if(F == 0) {
        *num_seams_out = 0;
        return NULL;
    }
    vector<vector<int>> dual_graph(F);
    unordered_map<long long, int> face_edge_map;
    for(int edge_idx = 0; edge_idx < topo->num_edges; edge_idx++) {
        int f0 = topo->edge_faces[2 * edge_idx + 0];
        int f1 = topo->edge_faces[2 * edge_idx + 1];

        if(f1 >= 0 && f0 >=0) {
            dual_graph[f0].push_back(f1);
            dual_graph[f1].push_back(f0);
            long long edge_key = make_edge_key(f0, f1);
            face_edge_map[edge_key] = edge_idx;
        } 
        if(f1 == -1 || f0 == -1) {
            boundry_edges_count++;
        }

    }
    // printf("Dual graph built with %d faces and %d boundary edges\n", F, boundry_edges_count);
    // STEP 2: Spanning tree via BFS
    set<int> tree_edges;
    vector<bool> visited(F, false);
    queue<int> q;
    // int components = 0;
    for(int start_face = 0; start_face < F; start_face++) {
        if(visited[start_face]) continue;
        // components++;
        q.push(start_face);
        visited[start_face] = true;

        while(!q.empty()) {
            int current_face = q.front();
            q.pop();

            for(int neighbor_face : dual_graph[current_face]) {
                if(!visited[neighbor_face]) {
                    visited[neighbor_face] = true;
                    q.push(neighbor_face);

                    long long edge_key = make_edge_key(current_face, neighbor_face);

                    if(face_edge_map.find(edge_key) != face_edge_map.end()) {
                        int edge_idx = face_edge_map[edge_key];
                        tree_edges.insert(edge_idx);
                    }
                }
            }
        }
    }

    // printf("Spanning tree completed with %d components\n", components);
    vector<float> vertex_defect(mesh->num_vertices, 0.0f);
    for(int v_idx = 0; v_idx < mesh->num_vertices; v_idx++) {
        vertex_defect[v_idx] = compute_angular_defect(mesh, v_idx);
    }

    // STEP 4: Initial seam candidates = non-tree interior edges with high angular defect
    seam_candidates.clear();
    vector<int> all_non_tree_edges; // Backup for smooth surfaces
    vector<pair<int, float>> defect_edges; // Store edges with their defects for open meshes
    
    // Adjust threshold based on whether mesh is open or closed
    // Open meshes (with boundaries) already have natural cuts, so be more conservative
    float threshold_multiplier = (boundry_edges_count > 0) ? 3.5f : 0.8f;
    
    for(int edge_idx = 0; edge_idx < topo->num_edges; edge_idx++) {
        int f0 = topo->edge_faces[2 * edge_idx + 0];
        int f1 = topo->edge_faces[2 * edge_idx + 1];

        if(f1 == -1 || f0 == -1) continue; // Skip boundary edges

        if(tree_edges.find(edge_idx) == tree_edges.end()) {
            // Not in tree = seam candidate
            // Check if either endpoint has significant angular defect
            all_non_tree_edges.push_back(edge_idx);
            
            int v0 = topo->edges[2 * edge_idx + 0];
            int v1 = topo->edges[2 * edge_idx + 1];
            
            float max_defect = max(fabsf(vertex_defect[v0]), fabsf(vertex_defect[v1]));
            
            // Only add if angular defect is significant
            if (max_defect > threshold_multiplier * angle_threshold_redians) {
                if (boundry_edges_count > 0) {
                    // For open meshes, collect edges with defects for ranking
                    defect_edges.push_back(make_pair(edge_idx, max_defect));
                } else {
                    // For closed meshes, add all high-defect edges
                    seam_candidates.insert(edge_idx);
                }
            }
        }
    }
    
    // For open meshes, limit to top 3 highest-defect edges
    if (boundry_edges_count > 0 && !defect_edges.empty()) {
        printf("  Open mesh detected - found %d high-defect edges:\n", (int)defect_edges.size());
        printf("  Threshold multiplier: %.2f, angle_threshold: %.4f rad (%.1f°)\n", 
               threshold_multiplier, angle_threshold_redians, angle_threshold);
        printf("  Effective threshold: %.4f rad\n", threshold_multiplier * angle_threshold_redians);
        
        sort(defect_edges.begin(), defect_edges.end(), 
             [](const pair<int, float>& a, const pair<int, float>& b) {
                 return a.second > b.second; // Sort by defect descending
             });
        
        for (size_t i = 0; i < defect_edges.size(); i++) {
            printf("    Edge %d: defect = %.4f rad (%.2f°) %s\n", 
                   defect_edges[i].first, 
                   defect_edges[i].second,
                   defect_edges[i].second * 180.0f / M_PI,
                   (i < 3) ? "[SELECTED]" : "[SKIPPED]");
        }
        
        int max_seams_open = min(3, (int)defect_edges.size());
        for (int i = 0; i < max_seams_open; i++) {
            seam_candidates.insert(defect_edges[i].first);
        }
    }

    // Fallback for smooth surfaces: add at least 1 seam to prevent severe distortion
    if (seam_candidates.empty() && !all_non_tree_edges.empty()) {
        // Add the first few non-tree edges as seams
        int num_fallback_seams = min((int)all_non_tree_edges.size(), max(1, F / 40));
        for (int i = 0; i < num_fallback_seams; i++) {
            seam_candidates.insert(all_non_tree_edges[i]);
        }
        printf("  (Added %d fallback seams for smooth surface)\n", num_fallback_seams);
    }

    printf("Total seam candidates selected: %zu\n", seam_candidates.size());



    // STEP 5: 
    // Convert to array
    // *num_seams_out = seam_candidates.size();
    
    *num_seams_out = static_cast<int>(seam_candidates.size());
    int* seams = (int*)malloc(*num_seams_out * sizeof(int));

    int idx = 0;
    for (int edge_idx : seam_candidates) {
        seams[idx++] = edge_idx;
    }

    printf("Detected %d seams\n", *num_seams_out);

    return seams;
}
