/* Copyright (c) 2018-2019 the `graphkernels` developers
 * All rights reserved.
 */

#ifndef GRAPHKERNELS_CPPKERNELS_REST_H_
#define GRAPHKERNELS_CPPKERNELS_REST_H_

#include <Eigen/Core>

#include <vector>

// Random Walk kernels
// ===================

Eigen::MatrixXd CalculateGeometricRandomWalkKernelPy(
        const std::vector<Eigen::MatrixXi>& E,
        const std::vector<std::vector<int>>& V_label,
        double par,
        int max_iterations,
        double eps);

Eigen::MatrixXd CalculateExponentialRandomWalkKernelPy(
        const std::vector<Eigen::MatrixXi>& E,
        const std::vector<std::vector<int>>& V_label,
        double par);

Eigen::MatrixXd CalculateKStepRandomWalkKernelPy(
        const std::vector<Eigen::MatrixXi>& E,
        const std::vector<std::vector<int>>& V_label,
        const std::vector<double>& par);

#endif  // GRAPHKERNELS_CPPKERNELS_REST_H_
