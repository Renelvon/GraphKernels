/* Copyright (c) 2018-2019 the `graphkernels` developers
 * All rigths reserved.
 */

#ifndef GRAPHKERNELS_CPPKERNELS_REST_H_
#define GRAPHKERNELS_CPPKERNELS_REST_H_

#include <Eigen/Core>

#include <vector>

using std::vector;

using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::MatrixXi;

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

Eigen::MatrixXd CalculateKernelPy(std::vector<Eigen::MatrixXi>& E,
                                  std::vector<std::vector<int>>& V_label,
                                  std::vector<double>& par,
                                  int kernel_type);

#endif  // GRAPHKERNELS_CPPKERNELS_REST_H_
