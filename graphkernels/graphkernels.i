/* Copyright (c) 2018-2019 the `graphkernels` developers
 * All rigths reserved.
 */

/* Define the SWIG module: graphkernels */

%module graphkernels

%{
/* line for specifying that the C file should be built as a python extension */
#define SWIG_FILE_WITH_INIT

#include <vector>

#include <Python.h>

#include "connected_graphlet.h"
#include "graphlet.h"
#include "rest.h"
#include "wl.h"
%}

// Include the built-in support for std::vector
%include <std_vector.i>
%include <numpy.i>
%include <eigen.i>


%init %{
import_array();
%}

namespace std {
%template(IntVector) std::vector<int>;
%template(IntIntVector) std::vector<std::vector<int>>;
%template(IntIntIntVector) std::vector<std::vector<std::vector<int>>>;
/*    %template(FloatVector) std::vector<float>; */
%template(DoubleVector) std::vector<double>;
%template(VecMatrixXi) std::vector<Eigen::MatrixXi>;
};

Eigen::MatrixXd WLKernelMatrix(std::vector<Eigen::MatrixXi>& E,
                               std::vector<std::vector<int>>& V_label,
                               std::vector<int>& num_v,
                               std::vector<int>& num_e,
                               std::vector<int>& degree_max,
                               int h_max);

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

Eigen::MatrixXd CalculateGraphletKernelPy(
    std::vector<std::vector<std::vector<int>>>& graph_adjlist_all,
    int k);

Eigen::MatrixXd CalculateConnectedGraphletKernelPy(
    std::vector<Eigen::MatrixXi>& graph_adj_all,
    std::vector<std::vector<std::vector<int>>>& graph_adjlist_all,
    int k);
