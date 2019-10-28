/* Copyright (c) 2018-2019 the `graphkernels` developers
 * All rigths reserved.
 */

#ifndef GRAPHKERNELS_CPPKERNELS_CONNECTED_GRAPHLET_H_
#define GRAPHKERNELS_CPPKERNELS_CONNECTED_GRAPHLET_H_

#include <Eigen/Core>

#include <vector>

// Connected Graphlet kernel for k = 3, 4, 5
// -----------------------------------------
Eigen::MatrixXd CalculateConnectedGraphletKernelPy(
    std::vector<Eigen::MatrixXi>& graph_adj_all,
    std::vector<std::vector<std::vector<int>>>& graph_adjlist_all,
    int k);

#endif  // GRAPHKERNELS_CPPKERNELS_CONNECTED_GRAPHLET_H_
