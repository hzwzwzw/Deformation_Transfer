/*
@Author: huang zewen

The project implements deformation transfer algorithm.

Included libraries: Eigen, OpenMesh
*/

#include "util.h"

/* Calculate v4 = v1 + (v2 − v1) × (v3 − v1)/√|(v2 − v1) × (v3 − v1)|*/
Point calculate_v4(const Point &v1, const Point &v2, const Point &v3)
{
    Point v4;
    Point v2_v1 = v2 - v1;
    Point v3_v1 = v3 - v1;
    Point cross = v2_v1.cross(v3_v1);
    double norm = cross.norm();
    v4 = v1 + cross / norm;
    return v4;
}

/* for each triangle face, calculate v4 and add it into mesh */
void add_v4(TriMesh &mesh)
{
    for (TriMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it)
    {
        TriMesh::FaceHandle fh = *f_it;
        TriMesh::FaceVertexIter fv_it = mesh.fv_iter(fh);
        Point v1 = mesh.point(*fv_it);
        Point v2 = mesh.point(*(++fv_it));
        Point v3 = mesh.point(*(++fv_it));
        Point v4 = calculate_v4(v1, v2, v3);
        mesh.add_vertex(v4);
    }
}

/* Main transformation function */
TriMesh deformation_transfer(TriMesh &s0, TriMesh &s1, TriMesh &t0)
{
}

/*
Input: s0, s1, t0
Output: t1
*/
int main(int argc, char *argv[])
{
    assert(argc == 5);
    std::string s0_path = argv[1];
    std::string s1_path = argv[2];
    std::string t0_path = argv[3];
    std::string t1_path = argv[4];

    TriMesh s0, s1, t0;
    OpenMesh::IO::read_mesh(s0, s0_path);
    OpenMesh::IO::read_mesh(s1, s1_path);
    OpenMesh::IO::read_mesh(t0, t0_path);

    TriMesh t1 = deformation_transfer(s0, s1, t0);

    // OpenMesh::IO::write_mesh(t1, t1_path);
    return 0;
}
