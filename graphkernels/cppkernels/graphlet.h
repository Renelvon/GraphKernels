/* Copyright (c) 2018-2019 the `graphkernels` developers
 * All rights reserved.
 */

#ifndef GRAPHKERNELS_CPPKERNELS_GRAPHLET_H_
#define GRAPHKERNELS_CPPKERNELS_GRAPHLET_H_

#include <Eigen/Core>

#include <vector>

// Graphlet kernel for k = 3, 4
// ----------------------------
Eigen::MatrixXd CalculateGraphletKernelPy(
        std::vector<std::vector<std::vector<int>>>& graph_adjlist_all,
        int k);

#endif  // GRAPHKERNELS_CPPKERNELS_GRAPHLET_H_
