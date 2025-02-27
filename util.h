#include <iostream>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <Eigen/Dense>

typedef OpenMesh::TriMesh_ArrayKernelT<> TriMesh;

void load_mesh(const char *filename, TriMesh &mesh)
{
    if (!OpenMesh::IO::read_mesh(mesh, filename))
    {
        std::cerr << "Error: Cannot read mesh file " << filename << std::endl;
        exit(1);
    }
}

void save_mesh(const char *filename, TriMesh &mesh)
{
    if (!OpenMesh::IO::write_mesh(mesh, filename))
    {
        std::cerr << "Error: Cannot write mesh file " << filename << std::endl;
        exit(1);
    }
}
