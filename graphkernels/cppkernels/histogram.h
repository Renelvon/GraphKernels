/* Copyright (c) 2018-2019 the `graphkernels` developers
 * All rigths reserved.
 */

#ifndef GRAPHKERNELS_CPPKERNELS_HISTOGRAM_H_
#define GRAPHKERNELS_CPPKERNELS_HISTOGRAM_H_

#include <Eigen/Core>

#include <vector>

double vertexHistogramKernel(std::vector<int>& v1_label,
                             std::vector<int>& v2_label,
                             double sigma);

double edgeHistogramKernel(Eigen::MatrixXi& e1,
                           Eigen::MatrixXi& e2,
                           double sigma);

double vertexEdgeHistogramKernel(Eigen::MatrixXi& e1,
                                 Eigen::MatrixXi& e2,
                                 std::vector<int>& v1_label,
                                 std::vector<int>& v2_label,
                                 double sigma);

double vertexVertexEdgeHistogramKernel(Eigen::MatrixXi& e1,
                                       Eigen::MatrixXi& e2,
                                       std::vector<int>& v1_label,
                                       std::vector<int>& v2_label,
                                       double lambda);

double computeKernelValue(Eigen::MatrixXi& e1,
                          Eigen::MatrixXi& e2,
                          std::vector<int>& v1_label,
                          std::vector<int>& v2_label,
                          std::vector<double>& par,
                          int kernel_type);

Eigen::MatrixXd CalculateKernelPy(std::vector<Eigen::MatrixXi>& E,
                                  std::vector<std::vector<int>>& V_label,
                                  std::vector<double>& par,
                                  int kernel_type);

#endif  // GRAPHKERNELS_CPPKERNELS_HISTOGRAM_H_
