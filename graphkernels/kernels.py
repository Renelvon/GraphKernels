"""
Functions for computing the graph kernels
"""

import numpy as np
from igraph import Graph

# FIXME: Avoid double-import by exporting names in __init__
from graphkernels import graphkernels as gkCpy

from .utilities import GetAdjMatList, GetGKInput


def _do_calculate_histogram(G, gk_par, kernel_id=None):
    if not isinstance(kernel_id, int):
        raise TypeError('Histogram Kernel index must be integer')

    if kernel_id <= 0:
        raise ValueError('Histogram Kernel index must be positive')

    # Extract graph info.
    E, V_label, _, _, _ = GetGKInput(G)

    # Compute designated kernel
    return gkCpy.CalculateHistogramKernelPy(E, V_label, gk_par, kernel_id)


# === Linear Kernels on Histograms ===


def CalculateEdgeHistKernel(G, par=-1.0):
    """Edge Histogram Kernel"""
    gk_par = gkCpy.DoubleVector([par])
    return _do_calculate_histogram(G, gk_par, 1)


def CalculateVertexHistKernel(G, par=-1.0):
    """Vertex Histogram Kernel"""
    gk_par = gkCpy.DoubleVector([par])
    return _do_calculate_histogram(G, gk_par, 2)


def CalculateVertexEdgeHistKernel(G, par=-1.0):
    """Vertex Edge Histogram Kernel"""
    gk_par = gkCpy.DoubleVector([par])
    return _do_calculate_histogram(G, gk_par, 3)


def CalculateVertexVertexEdgeHistKernel(G, par=1):
    """Vertex Vertex Edge Histogram Kernel"""
    gk_par = gkCpy.DoubleVector([par])
    return _do_calculate_histogram(G, gk_par, 4)


# === RBF Kernels on Histograms ===


def CalculateEdgeHistGaussKernel(G, par=1):
    """Edge Histogram RBF Kernel"""
    gk_par = gkCpy.DoubleVector([par])
    return _do_calculate_histogram(G, gk_par, 5)


def CalculateVertexHistGaussKernel(G, par=1):
    """Vertex Histogram RBF Kernel"""
    gk_par = gkCpy.DoubleVector([par])
    return _do_calculate_histogram(G, gk_par, 6)


def CalculateVertexEdgeHistGaussKernel(G, par=1):
    """Vertex Edge Histogram RBF Kernel"""
    gk_par = gkCpy.DoubleVector([par])
    return _do_calculate_histogram(G, gk_par, 7)


# === Random Walk Kernels ===


def CalculateGeometricRandomWalkKernel(
    G, par=1, max_iterations=100, eps=10.0 ** (-10)
):
    """Geometric Random Walk Kernel"""
    E, V_label, _, _, _ = GetGKInput(G)
    return gkCpy.CalculateGeometricRandomWalkKernelPy(
        E, V_label, par, max_iterations, eps
    )


def CalculateExponentialRandomWalkKernel(G, par=1):
    """Exponential Random Walk Kernel"""
    E, V_label, _, _, _ = GetGKInput(G)
    return gkCpy.CalculateExponentialRandomWalkKernelPy(E, V_label, par)


def CalculateKStepRandomWalkKernel(G, par=1):
    """K-step Random Walk Kernel"""
    # Allow user to provide own list of k-step weights
    if isinstance(par, (int, float, complex)):
        gk_par = gkCpy.DoubleVector([par])
    else:
        gk_par = gkCpy.DoubleVector(par)
    E, V_label, _, _, _ = GetGKInput(G)
    return gkCpy.CalculateKStepRandomWalkKernelPy(E, V_label, gk_par)


# === Advanced Kernels ===


def CalculateWLKernel(G, par=5):
    """Weisfeiler-Lehman Kernel

    Parametres
    ----------
    par : number of WL iterations
    """
    if not isinstance(par, int):
        raise TypeError('Number of WL iterations must be an integer')

    if par < 0:
        raise ValueError('Number of WL iterations must be non-negative')

    E, V_label, V_count, E_count, D_max = GetGKInput(G)  # Extract graph info.
    return gkCpy.WLKernelMatrix(E, V_label, V_count, E_count, D_max, par)


def CalculateGraphletKernel(G, par=4):
    """Graphlet Kernel

    Parametres
    ----------
    par : size of graphlets used (k)
    """
    if not isinstance(par, int):
        raise TypeError('Size of graphlets must be an integer')

    if par not in (3, 4):
        raise ValueError("Graphlet kernel supports only: k = 3 or 4")

    _, adj_list = GetAdjMatList(G)  # Extract graph info.
    return gkCpy.CalculateGraphletKernelPy(adj_list, par)


def CalculateConnectedGraphletKernel(G, par=4):
    """Connected Graphlet Kernel

    Parametres
    ----------
    par : size of graphlets used (k)
    """
    if not isinstance(par, int):
        raise TypeError('Size of graphlets must be an integer')

    if par not in (3, 4, 5):
        raise ValueError(
            "Connected Graphlet kernel supports only: k = 3, 4 or 5"
        )

    adj_mat, adj_list = GetAdjMatList(G)  # Extract graph info.
    return gkCpy.CalculateConnectedGraphletKernelPy(adj_mat, adj_list, par)


def _floyd_transform(gg):
    # TODO: Beautify.
    g_floyd_am = gg.shortest_paths_dijkstra()
    g_floyd_am = np.asarray(g_floyd_am).reshape(
        len(g_floyd_am), len(g_floyd_am)
    )
    g = Graph.Adjacency((g_floyd_am > 0).tolist())
    g.es['label'] = g_floyd_am[g_floyd_am.nonzero()]
    g.vs['id'] = np.arange(len(gg.vs['label']))
    g.vs['label'] = gg.vs['label']
    return g


def CalculateShortestPathKernel(G):
    """Shortest Path Kernel"""

    floyd_graphs = tuple(_floyd_transform(g) for g in G)
    G_floyd = np.array(floyd_graphs)

    return CalculateKStepRandomWalkKernel(G_floyd, par=(0, 1))
