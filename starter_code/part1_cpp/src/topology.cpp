/**
 * @file topology.cpp
 * @brief Topology builder implementation
 *
 * SKELETON - YOU IMPLEMENT THIS
 *
 * Algorithm:
 * 1. Extract all edges from triangles
 * 2. Ensure uniqueness (always store as v0 < v1)
 * 3. For each edge, find adjacent faces
 * 4. Validate using Euler characteristic
 */

#include "topology.h"
#include <stdlib.h>
#include <stdio.h>
#include <map>
#include <vector>
using namespace std;
/**
 * @brief Edge structure for uniqueness
 */
struct Edge {
    int v0, v1;

    Edge(int a, int b) {
        // Always store smaller vertex first
        if (a < b) {
            v0 = a;
            v1 = b;
        } else {
            v0 = b;
            v1 = a;
        }
    }

    bool operator<(const Edge& other) const {
        if (v0 != other.v0) return v0 < other.v0;
        return v1 < other.v1;
    }
};

/**
 * @brief Edge information
 */
struct EdgeInfo {
    int face0;
    int face1;

    EdgeInfo() : face0(-1), face1(-1) {}
};

TopologyInfo* build_topology(const Mesh* mesh) {
    if (!mesh) return NULL;

    map<Edge, EdgeInfo> edge_map;

    for(int face_idx = 0; face_idx < mesh->num_triangles; face_idx++) {
        int v0 = mesh->triangles[3 * face_idx + 0];
        int v1 = mesh->triangles[3 * face_idx + 1];
        int v2 = mesh->triangles[3 * face_idx + 2];

        Edge edges[3] = { Edge(v0, v1), Edge(v1, v2), Edge(v2, v0) };

        for(int e = 0; e < 3; e++) {
            Edge& edge = edges[e];
            
            if(edge_map.find(edge) == edge_map.end()) {
                edge_map[edge].face0 = face_idx;
            } else {
                edge_map[edge].face1 = face_idx;
            }
        }
    }

    // TODO: Implement topology building
    //
    // Steps:
    // 1. Create std::map<Edge, EdgeInfo> to collect edges
    // 2. Iterate through all triangles
    //    For each triangle, extract 3 edges
    //    Add to map, tracking which faces use each edge
    // 3. Convert map to arrays (edges, edge_faces)
    // 4. Allocate TopologyInfo and fill arrays
    //
    // Hints:
    // - Use Edge struct for automatic ordering
    // - Each edge should have 1 or 2 adjacent faces
    // - Boundary edges have only 1 face (set face1 = -1)
    //
    // See reference/topology_example.cpp for complete example



    TopologyInfo* topo = (TopologyInfo*)malloc(sizeof(TopologyInfo));

    // Initialize to safe defaults (prevents crashes before implementation)
    topo->edges = NULL;
    topo->num_edges = 0;
    topo->edge_faces = NULL;

    // TODO: Your implementation here
    int num_edges = edge_map.size();
    topo->num_edges = num_edges;

    topo->edges = (int*)malloc(sizeof(int) * 2 * num_edges);
    topo->edge_faces = (int*)malloc(sizeof(int) * 2 * num_edges);

    int edge_idx = 0;
    for(auto& pair : edge_map) {
        Edge edge = pair.first;
        EdgeInfo info = pair.second;

        topo->edges[2 * edge_idx + 0] = edge.v0;
        topo->edges[2 * edge_idx + 1] = edge.v1;

        topo->edge_faces[2 * edge_idx + 0] = info.face0;
        topo->edge_faces[2 * edge_idx + 1] = info.face1;

        edge_idx++;
    }

    validate_topology(mesh, topo);

    return topo;
}

void free_topology(TopologyInfo* topo) {
    if (!topo) return;

    if (topo->edges) free(topo->edges);
    if (topo->edge_faces) free(topo->edge_faces);
    free(topo);
}

int validate_topology(const Mesh* mesh, const TopologyInfo* topo) {
    if (!mesh || !topo) return 0;

    int V = mesh->num_vertices;
    int E = topo->num_edges;
    int F = mesh->num_triangles;

    int euler = V - E + F;

    printf("Topology validation:\n");
    printf("  V=%d, E=%d, F=%d\n", V, E, F);
    printf("  Euler characteristic: %d (expected 2 for closed mesh)\n", euler);

    // Closed meshes should have Euler = 2
    // Open meshes or meshes with holes may differ
    if (euler != 2) {
        printf("  Warning: Non-standard Euler characteristic\n");
        printf("  (This may be OK for open meshes or meshes with boundaries)\n");
    }

    return 1;
}
