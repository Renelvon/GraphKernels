/* Copyright (c) 2018-2019 the `graphkernels` developers
 * All rigths reserved.
 */

/* Define the SWIG module: graphkernels */

%module graphkernels

%{
/* line for specifying that the C file should be built as a python extension */
#include <memory>
#include <vector>
#define SWIG_FILE_WITH_INIT
#include <Python.h>
#define ITPP_EXPORT
#include "graphkernels.h"
%}

// Include the built-in support for std::vector
%include <typemaps.i>
%include <std_vector.i>
%include <std_string.i>
%include <std_shared_ptr.i>
%include <numpy.i>
%include <eigen.i>
%include <stl.i>



%init %{
import_array();
%}

using namespace std;

namespace std {
%template(IntVector) vector<int>;
%template(IntIntVector) vector<vector<int>>;
%template(IntIntIntVector) vector<vector<vector<int>>>;
/*    %template(FloatVector) vector<float>; */
%template(DoubleVector) vector<double>;
%template(VecMatrixXi) vector<Eigen::MatrixXi>;
};

// All the kernel functions
double edgeHistogramKernel(Eigen::MatrixXi& e1,
                           Eigen::MatrixXi& e2,
                           double sigma);
double vertexHistogramKernel(vector<int>& v1_label,
                             vector<int>& v2_label,
                             double sigma);
double vertexEdgeHistogramKernel(Eigen::MatrixXi& e1,
                                 Eigen::MatrixXi& e2,
                                 vector<int>& v1_label,
                                 vector<int>& v2_label,
                                 double sigma);
double vertexVertexEdgeHistogramKernel(Eigen::MatrixXi& e1,
                                       Eigen::MatrixXi& e2,
                                       vector<int>& v1_label,
                                       vector<int>& v2_label,
                                       double lambda);

double geometricRandomWalkKernel(Eigen::MatrixXi& e1,
                                 Eigen::MatrixXi& e2,
                                 vector<int>& v1_label,
                                 vector<int>& v2_label,
                                 double lambda);
double exponentialRandomWalkKernel(Eigen::MatrixXi& e1,
                                   Eigen::MatrixXi& e2,
                                   vector<int>& v1_label,
                                   vector<int>& v2_label,
                                   double beta);
double kstepRandomWalkKernel(Eigen::MatrixXi& e1,
                             Eigen::MatrixXi& e2,
                             vector<int>& v1_label,
                             vector<int>& v2_label,
                             vector<double>& lambda_list);

Eigen::MatrixXd WLKernelMatrix(vector<Eigen::MatrixXi>& E,
                               vector<vector<int>>& V_label,
                               vector<int>& num_v,
                               vector<int>& num_e,
                               vector<int>& degree_max,
                               int h_max);

double computeKernelValue(Eigen::MatrixXi& e1,
                          Eigen::MatrixXi& e2,
                          vector<int>& v1_label,
                          vector<int>& v2_label,
                          vector<double>& par,
                          int kernel_type);

Eigen::MatrixXd CalculateKernelPy(vector<Eigen::MatrixXi>& E,
                                  vector<vector<int>>& V_label,
                                  vector<double>& par,
                                  int kernel_type);

Eigen::VectorXd countGraphletsFour(vector<vector<int>>& al,
                                   Eigen::VectorXd& count_gr);

Eigen::VectorXd countGraphletsThree(vector<vector<int>>& al,
                                    Eigen::VectorXd& count_gr);

Eigen::VectorXd countConnectedGraphletsFive(Eigen::MatrixXi& am,
                                            vector<vector<int>>& al,
                                            Eigen::VectorXd& count_gr);

Eigen::VectorXd countConnectedGraphletsFour(Eigen::MatrixXi& am,
                                            vector<vector<int>>& al,
                                            Eigen::VectorXd& count_gr);

Eigen::VectorXd countConnectedGraphletsThree(Eigen::MatrixXi& am,
                                             vector<vector<int>>& al,
                                             Eigen::VectorXd& count_gr);

Eigen::MatrixXd CalculateGraphletKernelPy(
    vector<vector<vector<int>>>& graph_adjlist_all,
    int k);

Eigen::MatrixXd CalculateConnectedGraphletKernelPy(
    vector<Eigen::MatrixXi>& graph_adj_all,
    vector<vector<vector<int>>>& graph_adjlist_all,
    int k);
