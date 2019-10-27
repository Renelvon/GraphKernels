/* Copyright (c) 2018-2019 the `graphkernels` developers
 * All rigths reserved.
 */

#ifndef GRAPHKERNELS_CPPKERNELS_GRAPHLET_H_
#define GRAPHKERNELS_CPPKERNELS_GRAPHLET_H_

#include <Eigen/Core>

#include <vector>

using std::vector;

using Eigen::MatrixXd;

// Graphlet kernel for k = 3, 4
// ----------------------------
MatrixXd CalculateGraphletKernelPy(
    vector<vector<vector<int>>>& graph_adjlist_all,
    int k);

#endif  // GRAPHKERNELS_CPPKERNELS_GRAPHLET_H_
