/* Copyright (c) 2018-2019 the `graphkernels` developers
 * All rigths reserved.
 */

#ifndef GKEXTCPY_GKEXTCPY_H_
#define GKEXTCPY_GKEXTCPY_H_

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

double selectLinearGaussian(vector<int>& h1, vector<int>& h2, double sigma);

int productMapping(MatrixXi& e1,
                   MatrixXi& e2,
                   vector<int>& v1_label,
                   vector<int>& v2_label,
                   MatrixXi& H);

MatrixXd productAdjacency(MatrixXi& e1,
                          MatrixXi& e2,
                          vector<int>& v1_label,
                          vector<int>& v2_label,
                          MatrixXi& H);

void bucketsort(vector<int>& x, vector<int>& index, int label_max);

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
                           vector<int>& V_count,
                           vector<int>& E_count,
                           vector<int>& D_max,
                           vector<double>& par,
                           int kernel_type);

// Graphlet kernels
// ================

// Auxiliary functions
// -------------------
int find_min(int a, int b, int c);

void card_ThreeInter(vector<int>& L1,
                     vector<int>& L2,
                     vector<int>& L3,
                     vector<int>& card);

void getIndices(vector<int>& o_set1,
                vector<int>& o_set2,
                vector<int>& inter,
                vector<int>& diff1,
                vector<int>& diff2);


void getCardinality(vector<int>& o_set1,
                    vector<int>& o_set2,
                    vector<double>& card);

void getMinValue(MatrixXi& iam, vector<int>& idx, vector<int>& sums);

// Graphlet kernel for k = 3, 4
// ----------------------------
VectorXd countGraphletsThree(vector<vector<int>>& al, VectorXd& count_gr);

VectorXd countGraphletsFour(vector<vector<int>>& al, VectorXd& count_gr);

MatrixXd CalculateGraphletKernelPy(
    vector<MatrixXi>& graph_adj_all,
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

#endif  // GKEXTCPY_GKEXTCPY_H_
