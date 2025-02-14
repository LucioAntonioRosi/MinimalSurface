#include <fstream>
#include <iostream>
#include <sys/stat.h>
#include <sys/types.h>
#include "export_mesh.h"

// Function to create directory if it does not exist
void create_directory(const std::string& directory) {
    struct stat info;
    if (stat(directory.c_str(), &info) != 0) {
        // Directory does not exist, create it
        if (mkdir(directory.c_str(), 0777) == -1) {
            std::cerr << "Error: Cannot create directory " << directory << "\n";
        }
    } else if (!(info.st_mode & S_IFDIR)) {
        std::cerr << "Error: " << directory << " is not a directory!\n";
    }
}

void export_mesh_to_csv(const Mesh& mesh, const TArray<double>& solution, const std::string& filename) {
    // Create postProcessing directory if it does not exist
    create_directory("postProcessing");

    std::ofstream file("postProcessing/" + filename);
    
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file for writing!\n";
        return;
    }
    
    // Write header
    file << "x,y,z,type\n";

    // Save vertices with solution values
    for (int i = 0; i < mesh.vtx_count; ++i) {
        file << mesh.vertices[i].x << "," << mesh.vertices[i].y << "," << solution[i] << ",vertex\n";
    }

    // Save triangles (edges for visualization)
    for (int i = 0; i < mesh.tri_count; ++i) {
        int v0 = mesh.triangles[i].x;
        int v1 = mesh.triangles[i].y;
        int v2 = mesh.triangles[i].z;
        file << mesh.vertices[v0].x << "," << mesh.vertices[v0].y << "," << solution[v0] << ",edge\n";
        file << mesh.vertices[v1].x << "," << mesh.vertices[v1].y << "," << solution[v1] << ",edge\n";
        file << "\n";
        file << mesh.vertices[v1].x << "," << mesh.vertices[v1].y << "," << solution[v1] << ",edge\n";
        file << mesh.vertices[v2].x << "," << mesh.vertices[v2].y << "," << solution[v2] << ",edge\n";
        file << "\n";
        file << mesh.vertices[v2].x << "," << mesh.vertices[v2].y << "," << solution[v2] << ",edge\n";
        file << mesh.vertices[v0].x << "," << mesh.vertices[v0].y << "," << solution[v0] << ",edge\n";
        file << "\n";
    }

    file.close();
}

void export_mesh_to_csv(const Mesh& mesh, const TVec3<TArray<double>>& solution, const std::string& filename) {
    // Create postProcessing directory if it does not exist
    create_directory("postProcessing");

    std::ofstream file("postProcessing/" + filename);
    
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file for writing!\n";
        return;
    }
    
    // Write header
    file << "u1,u2,u3,type\n";

    // Save vertices with solution values
    for (int i = 0; i < mesh.vtx_count; ++i) {
        file << solution[0][i] << "," << solution[1][i] << "," << solution[2][i] << ",vertex\n";
    }

    // Save triangles (edges for visualization)
    for (int i = 0; i < mesh.tri_count; ++i) {
        int v0 = mesh.triangles[i].x;
        int v1 = mesh.triangles[i].y;
        int v2 = mesh.triangles[i].z;
        file << solution[0][v0] << "," << solution[1][v0] << "," << solution[2][v0] << ",edge\n";
        file << solution[0][v1] << "," << solution[1][v1] << "," << solution[2][v1] << ",edge\n";
        file << "\n";
        file << solution[0][v1] << "," << solution[1][v1] << "," << solution[2][v1] << ",edge\n";
        file << solution[0][v2] << "," << solution[1][v2] << "," << solution[2][v2] << ",edge\n";
        file << "\n";
        file << solution[0][v2] << "," << solution[1][v2] << "," << solution[2][v2] << ",edge\n";
        file << solution[0][v0] << "," << solution[1][v0] << "," << solution[2][v0] << ",edge\n";
        file << "\n";
    }

    file.close();
}