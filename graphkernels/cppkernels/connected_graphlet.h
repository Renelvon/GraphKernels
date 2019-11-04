/* Copyright (c) 2018-2019 the `graphkernels` developers
 * All rights reserved.
 */

#ifndef GRAPHKERNELS_CPPKERNELS_CONNECTED_GRAPHLET_H_
#define GRAPHKERNELS_CPPKERNELS_CONNECTED_GRAPHLET_H_

#include <Eigen/Core>

#include <vector>

Eigen::MatrixXd CalculateConnectedGraphletKernelThreePy(
        const std::vector<Eigen::MatrixXi>& graph_adj_all,
        const std::vector<std::vector<std::vector<int>>>& graph_adjlist_all);

Eigen::MatrixXd CalculateConnectedGraphletKernelFourPy(
        const std::vector<Eigen::MatrixXi>& graph_adj_all,
        const std::vector<std::vector<std::vector<int>>>& graph_adjlist_all);

Eigen::MatrixXd CalculateConnectedGraphletKernelFivePy(
        const std::vector<Eigen::MatrixXi>& graph_adj_all,
        const std::vector<std::vector<std::vector<int>>>& graph_adjlist_all);

#endif  // GRAPHKERNELS_CPPKERNELS_CONNECTED_GRAPHLET_H_
