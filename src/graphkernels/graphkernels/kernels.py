"""
Functions for computing the graph kernels
"""

import numpy as np
from igraph import Graph

import GKextCPy as gkCpy

from .utilities import GetGKInput, GetAdjMatList


def _do_calculate(G, gk_par, kernel_id=None):
    if not isinstance(kernel_id, int):
        raise TypeError('Kernel index must be integer')

    if kernel_id <= 0:
        raise ValueError('Kernel index must be positive')

    # Extract graph info.
    E, V_label, V_count, E_count, D_max = GetGKInput(G)

    # Compute designated kernel
    return gkCpy.CalculateKernelPy(
        E, V_label, V_count, E_count, D_max,
        gk_par,
        kernel_id
    )


# === Linear Kernels on Histograms ===

def CalculateEdgeHistKernel(G, par=-1.0):
    """Edge Histogram Kernel"""
    gk_par = gkCpy.DoubleVector([par])
    return _do_calculate(G, gk_par, 1)


def CalculateVertexHistKernel(G, par=-1.0):
    """Vertex Histogram Kernel"""
    gk_par = gkCpy.DoubleVector([par])
    return _do_calculate(G, gk_par, 2)


def CalculateVertexEdgeHistKernel(G, par=-1.0):
    """Vertex Edge Histogram Kernel"""
    gk_par = gkCpy.DoubleVector([par])
    return _do_calculate(G, gk_par, 3)


def CalculateVertexVertexEdgeHistKernel(G, par=1):
    """Vertex Vertex Edge Histogram Kernel"""
    gk_par = gkCpy.DoubleVector([par])
    return _do_calculate(G, gk_par, 4)


# === RBF Kernels on Histograms ===

def CalculateEdgeHistGaussKernel(G, par=1):
    """Edge Histogram RBF Kernel"""
    gk_par = gkCpy.DoubleVector([par])
    return _do_calculate(G, gk_par, 5)


def CalculateVertexHistGaussKernel(G, par=1):
    """Vertex Histogram RBF Kernel"""
    gk_par = gkCpy.DoubleVector([par])
    return _do_calculate(G, gk_par, 6)


def CalculateVertexEdgeHistGaussKernel(G, par=1):
    """Vertex Edge Histogram RBF Kernel"""
    gk_par = gkCpy.DoubleVector([par])
    return _do_calculate(G, gk_par, 7)


# === Random Walk Kernels ===

def CalculateGeometricRandomWalkKernel(G, par=1):
    """Geometric Random Walk Kernel"""
    gk_par = gkCpy.DoubleVector([par])
    return _do_calculate(G, gk_par, 8)


def CalculateExponentialRandomWalkKernel(G, par=1):
    """Exponential Random Walk Kernel"""
    gk_par = gkCpy.DoubleVector([par])
    return _do_calculate(G, gk_par, 9)


def CalculateKStepRandomWalkKernel(G, par=1):
    """K-step Random Walk Kernel"""
    # Allow user to provide own list of k-step weights
    if isinstance(par, (int, float, complex)):
        par = [par]
    gk_par = gkCpy.DoubleVector(par)
    return _do_calculate(G, gk_par, 10)


# === Advanced Kernels ===

def CalculateWLKernel(G, par=5):
    """Weisfeiler-Lehman Kernel Kernel"""

    # Extract graph info
    E, V_label, V_count, E_count, D_max = GetGKInput(G)

    #par = nuber of WL iterations
    par = int(par)

    return gkCpy.WLKernelMatrix(E, V_label, V_count, E_count, D_max, par)


def CalculateGraphletKernel(G, par=4):
    """Graphlet Kernel"""

    # If k<3 then assign k=3
    if par < 3:

        par = 3
        print("Warning: k=3 is used (k = 3 or 4 is supported)")

    # If k>4 then assign k=4
    if par > 4:

        par = 4
        print("Warning: k=4 is used (k = 3 or 4 is supported)")

    # Extract graph info
    adj_mat, adj_list = GetAdjMatList(G)
    par = int(par)

    return gkCpy.CalculateGraphletKernelPy(adj_mat, adj_list, par)


def CalculateConnectedGraphletKernel(G, par=4):
    """Connected Graphlet Kernel"""

    # If k<3 then assign k=3
    if par < 3:

        par = 3
        print("Warning: k=3 is used (k = 3, 4 or 5 is supported)")

    # If k>5 then assign k=5
    if par > 5:

        par = 5
        print("Warning: k=5 is used (k = 3, 4 or 5 is supported)")

    # Extract graph info
    adj_mat, adj_list = GetAdjMatList(G)

    return gkCpy.CalculateConnectedGraphletKernelPy(adj_mat, adj_list, par)


def CalculateShortestPathKernel(G):
    """Shortest Path Kernel"""

    G_floyd = []
    for i in range(len(G)):

        g_floyd_am = G[i].shortest_paths_dijkstra()
        g_floyd_am = np.asarray(g_floyd_am).reshape(len(g_floyd_am), len(g_floyd_am))
        g = Graph.Adjacency((g_floyd_am > 0).tolist())
        g.es['label'] = g_floyd_am[g_floyd_am.nonzero()]
        g.vs['id'] = np.arange(len(G[i].vs['label']))
        g.vs['label'] = G[i].vs['label']
        G_floyd.append(g)

    G_floyd = np.array(G_floyd)

    return CalculateKStepRandomWalkKernel(G_floyd, par=(0, 1))
