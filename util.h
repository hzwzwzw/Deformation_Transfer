#include <iostream>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <Eigen/Dense>
#include <vector>

typedef OpenMesh::TriMesh_ArrayKernelT<> TriMesh;
typedef OpenMesh::DefaultTraits::Point Point;

void LUdecomposition(Eigen::MatrixXd &A, Eigen::MatrixXd &L, Eigen::MatrixXd &U)
{
    int n = A.rows();
    L = Eigen::MatrixXd::Zero(n, n);
    U = Eigen::MatrixXd::Zero(n, n);
    for (int i = 0; i < n; i++)
    {
        L(i, i) = 1;
    }
    for (int i = 0; i < n; i++)
    {
        for (int j = i; j < n; j++)
        {
            U(i, j) = A(i, j);
            for (int k = 0; k < i; k++)
            {
                U(i, j) -= L(i, k) * U(k, j);
            }
        }
        for (int j = i + 1; j < n; j++)
        {
            L(j, i) = A(j, i);
            for (int k = 0; k < i; k++)
            {
                L(j, i) -= L(j, k) * U(k, i);
            }
            L(j, i) /= U(i, i);
        }
    }
}

void LUsolve(Eigen::MatrixXd &L, Eigen::MatrixXd &U, Eigen::VectorXd &b, Eigen::VectorXd &x)
{
    int n = L.rows();
    Eigen::VectorXd y = Eigen::VectorXd::Zero(n);
    for (int i = 0; i < n; i++)
    {
        y(i) = b(i);
        for (int j = 0; j < i; j++)
        {
            y(i) -= L(i, j) * y(j);
        }
    }
    for (int i = n - 1; i >= 0; i--)
    {
        x(i) = y(i);
        for (int j = i + 1; j < n; j++)
        {
            x(i) -= U(i, j) * x(j);
        }
        x(i) /= U(i, i);
    }
}
