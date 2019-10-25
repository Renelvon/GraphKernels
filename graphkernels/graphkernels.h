/* Copyright (c) 2018-2019 the `graphkernels` developers
 * All rigths reserved.
 */

#ifndef GRAPHKERNELS_GRAPHKERNELS_H_
#define GRAPHKERNELS_GRAPHKERNELS_H_

#include <algorithm>
#include <functional>
#include <numeric>
#include <set>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <Eigen/LU>
#include <Eigen/Sparse>

#define Int int32_t

typedef Eigen::Triplet<double> T;


using std::accumulate;
using std::greater;
using std::iota;
using std::max;
using std::min;
using std::set;
using std::vector;

using Eigen::FullPivLU;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::SelfAdjointEigenSolver;
using Eigen::SparseMatrix;

// Simple Kernels
// ==============

double vertexHistogramKernel(vector<int>& v1_label,
                             vector<int>& v2_label,
                             double sigma);

double edgeHistogramKernel(MatrixXi& e1, MatrixXi& e2, double sigma);

double vertexEdgeHistogramKernel(MatrixXi& e1,
                                 MatrixXi& e2,
                                 vector<int>& v1_label,
                                 vector<int>& v2_label,
                                 double sigma);

double vertexVertexEdgeHistogramKernel(MatrixXi& e1,
                                       MatrixXi& e2,
                                       vector<int>& v1_label,
                                       vector<int>& v2_label,
                                       double lambda);

// Random Walk kernels
// ===================

double geometricRandomWalkKernel(MatrixXi& e1,
                                 MatrixXi& e2,
                                 vector<int>& v1_label,
                                 vector<int>& v2_label,
                                 double lambda);

double exponentialRandomWalkKernel(MatrixXi& e1,
                                   MatrixXi& e2,
                                   vector<int>& v1_label,
                                   vector<int>& v2_label,
                                   double beta);

double kstepRandomWalkKernel(MatrixXi& e1,
                             MatrixXi& e2,
                             vector<int>& v1_label,
                             vector<int>& v2_label,
                             vector<double>& lambda_list);

double computeKernelValue(MatrixXi& e1,
                          MatrixXi& e2,
                          vector<int>& v1_label,
                          vector<int>& v2_label,
                          vector<double>& par,
                          int kernel_type);

// WL kernel
// =========

MatrixXd WLKernelMatrix(vector<MatrixXi>& E,
                        vector<vector<int>>& V_label,
                        vector<int>& num_v,
                        vector<int>& num_e,
                        vector<int>& degree_max,
                        int h_max);

MatrixXd CalculateKernelPy(vector<MatrixXi>& E,
                           vector<vector<int>>& V_label,
                           vector<double>& par,
                           int kernel_type);

// Graphlet kernels
// ================

// Graphlet kernel for k = 3, 4
// ----------------------------
VectorXd countGraphletsThree(vector<vector<int>>& al, VectorXd& count_gr);

VectorXd countGraphletsFour(vector<vector<int>>& al, VectorXd& count_gr);

MatrixXd CalculateGraphletKernelPy(
    vector<vector<vector<int>>>& graph_adjlist_all,
    int k);

// Connected Graphlet kernel for k = 3, 4, 5
// -----------------------------------------

VectorXd countConnectedGraphletsThree(MatrixXi& am,
                                      vector<vector<int>>& al,
                                      VectorXd& count_gr);

VectorXd countConnectedGraphletsFour(MatrixXi& am,
                                     vector<vector<int>>& al,
                                     VectorXd& count_gr);

VectorXd countConnectedGraphletsFive(MatrixXi& am,
                                     vector<vector<int>>& al,
                                     VectorXd& count_gr);

MatrixXd CalculateConnectedGraphletKernelPy(
    vector<MatrixXi>& graph_adj_all,
    vector<vector<vector<int>>>& graph_adjlist_all,
    int k);

#endif  // GRAPHKERNELS_GRAPHKERNELS_H_
