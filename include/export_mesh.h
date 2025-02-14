#pragma once

#include "my_mesh.h"
#include "vec3.h"
#include "array.h"
#include <fstream>

void export_mesh_to_csv(const Mesh& mesh, const TArray<double>& solution, const std::string& filename);
void export_mesh_to_csv(const Mesh& mesh, const TVec3<TArray<double>>& solution, const std::string& filename);
