/* Copyright (c) 2018-2019 the `graphkernels` developers
 * All rigths reserved.
 */

#ifndef GRAPHKERNELS_CPPKERNELS_WL_H_
#define GRAPHKERNELS_CPPKERNELS_WL_H_

#include <Eigen/Core>

#include <vector>

using std::vector;

using Eigen::MatrixXd;
using Eigen::MatrixXi;

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

#endif  // GRAPHKERNELS_CPPKERNELS_WL_H_
