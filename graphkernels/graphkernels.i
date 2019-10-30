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
#include "histogram.h"
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

// Instantiate vector templates.
namespace std {
    %template(IntVector) std::vector<int>;
    %template(IntIntVector) std::vector<std::vector<int>>;
    %template(IntIntIntVector) std::vector<std::vector<std::vector<int>>>;
    // %template(FloatVector) std::vector<float>;
    %template(DoubleVector) std::vector<double>;
    %template(VecMatrixXi) std::vector<Eigen::MatrixXi>;
};

// Include header files with kernel prototypes
%include "connected_graphlet.h"
%include "graphlet.h"
%include "histogram.h"
%include "rest.h"
%include "wl.h"
