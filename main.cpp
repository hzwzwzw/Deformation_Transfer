/*
@Author: huang zewen

The project implements deformation transfer algorithm.

Included libraries: Eigen, OpenMesh
*/

#include "util.h"

int get_vertex_index(TriMesh::VertexHandle vh)
{
    return vh.idx();
}

int get_face_index(TriMesh::FaceHandle fh)
{
    return fh.idx();
}

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

/* Calculate Vhat_i_inv */
std::vector<Eigen::MatrixXd> calculate_matrixVhati_inv(TriMesh &mesh)
{
    int newindex = mesh.n_vertices();
    const int n = mesh.n_vertices();
    const int m = mesh.n_faces();
    std::vector<Eigen::MatrixXd> Vhati_inv;
    std::cout << "n:" << n << " m:" << m << std::endl;
    for (TriMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it)
    {
        TriMesh::FaceHandle fh = *f_it;
        TriMesh::FaceVertexIter fv_it = mesh.fv_iter(fh);
        Point v1 = mesh.point(*fv_it);
        int index1 = get_vertex_index(*fv_it);
        Point v2 = mesh.point(*(++fv_it));
        int index2 = get_vertex_index(*fv_it);
        Point v3 = mesh.point(*(++fv_it));
        int index3 = get_vertex_index(*fv_it);
        Point v4 = calculate_v4(v1, v2, v3);
        int index4 = newindex++; // or n + get_face_index(fh)
        // create matrix (n + m) * 3
        // where the row index1 is -1, -1, -1
        // the row index2 is 1, 0, 0
        // the row index3 is 0, 1, 0
        // the row index4 is 0, 0, 1
        Eigen::MatrixXd A(n + m, 3);
        A.setZero();
        // std::cout << "index:" << index1 << " " << index2 << " " << index3 << " " << index4 << std::endl;
        A(index1, 0) = -1;
        A(index1, 1) = -1;
        A(index1, 2) = -1;
        A(index2, 0) = 1;
        A(index3, 1) = 1;
        A(index4, 2) = 1;
        Eigen::Matrix3d m;
        m << v2[0] - v1[0], v3[0] - v1[0], v4[0] - v1[0],
            v2[1] - v1[1], v3[1] - v1[1], v4[1] - v1[1],
            v2[2] - v1[2], v3[2] - v1[2], v4[2] - v1[2];
        // A * m^-1
        Vhati_inv.push_back(A * m.inverse());
    }
    std::cout << "Vhati_inv size: " << Vhati_inv.size() << std::endl;
    return Vhati_inv;
}

/* Main transformation function */
TriMesh deformation_transfer(TriMesh &s0, TriMesh &s1, TriMesh &t0)
{
    // here we assume s0, s1, t0 have the same number of vertices and faces, or same topology
    // so we can use the same face handle to iterate over the three meshes
    // and we only need to implement the vertex formulation of deformation transfer, that is:
    // E = ∑_{i=1}^{M}||S_{sj}-T_{tj}||^2 = ||A_Q X^T - B||^2, where
    // AQ∈R3|M|×(n+m), bQ∈R3|M|×3.
    // the final result is: ATAx=ATb.

    const int n = s0.n_vertices();
    const int m = s0.n_faces();

    // compute Vhat_i^-1 , shape: each (m + n) * 3
    std::vector<Eigen::MatrixXd> Vhati_inv;
    Vhati_inv = calculate_matrixVhati_inv(s0);

    // compute target function of deformation transfer
    // A_Q = -(Vhat_1 ^ -1) ^ T\\ (Vhat_2 ^ -1) ^ T\\ ...\\ (Vhat_M ^ -1) ^ T
    Eigen::MatrixXd A_Q(3 * m, n + m);
    A_Q.setZero();
    const int width = n + m;
    for (int i = 0; i < s0.n_faces(); i++)
    {
        Eigen::MatrixXd Vhati = -Vhati_inv[i].transpose();
        A_Q.block(3 * i, 0, 3, width) = Vhati;
    }
    // std::cout << A_Q;
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

    std::cout << "s0: " << s0.n_vertices() << " vertices, " << s0.n_faces() << " faces" << std::endl;
    std::cout << "s1: " << s1.n_vertices() << " vertices, " << s1.n_faces() << " faces" << std::endl;
    std::cout << "t0: " << t0.n_vertices() << " vertices, " << t0.n_faces() << " faces" << std::endl;

    TriMesh t1 = deformation_transfer(s0, s1, t0);

    // OpenMesh::IO::write_mesh(t1, t1_path);
    return 0;
}
