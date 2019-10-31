"""
Functions for computing the graph kernels
"""

import collections
import warnings

import numpy as np
from igraph import Graph

# FIXME: Avoid double-import by exporting names in __init__
from graphkernels import graphkernels as gkCpy

from .utilities import GetAdjMatList, GetGKInput


# === Linear Kernels on Histograms ===


def CalculateEdgeHistKernel(G):
    """Edge Histogram Kernel"""
    E, V_label, _, _, _ = GetGKInput(G)
    return gkCpy.CalculateHistogramKernelPy(E, V_label, -1.0, 1)


def CalculateVertexHistKernel(G):
    """Vertex Histogram Kernel"""
    E, V_label, _, _, _ = GetGKInput(G)
    return gkCpy.CalculateHistogramKernelPy(E, V_label, -1.0, 2)


def CalculateVertexEdgeHistKernel(G):
    """Vertex Edge Histogram Kernel"""
    E, V_label, _, _, _ = GetGKInput(G)
    return gkCpy.CalculateHistogramKernelPy(E, V_label, -1.0, 3)


def CalculateVertexVertexEdgeHistKernel(G, par=1.0):
    """Vertex Vertex Edge Histogram Kernel"""
    if not isinstance(par, (float, int)):
        raise TypeError('par must be a scalar (float or integer)')

    if par == 0:
        warnings.warn('Invoking kernel with par == 0.0')

    E, V_label, _, _, _ = GetGKInput(G)
    return gkCpy.CalculateHistogramKernelPy(E, V_label, float(par), 4)


# === RBF Kernels on Histograms ===


def CalculateEdgeHistGaussKernel(G, gamma=0.5):
    """Edge Histogram RBF Kernel"""
    if not isinstance(gamma, (float, int)):
        raise TypeError('gamma must be a positive scalar (float or integer)')

    if gamma <= 0.0:
        raise ValueError('gamma must be a positive scalar (float or integer)')

    E, V_label, _, _, _ = GetGKInput(G)
    return gkCpy.CalculateHistogramKernelPy(E, V_label, float(gamma), 5)


def CalculateVertexHistGaussKernel(G, gamma=0.5):
    """Vertex Histogram RBF Kernel"""
    if not isinstance(gamma, (float, int)):
        raise TypeError('gamma must be a positive scalar (float or integer)')

    if gamma <= 0.0:
        raise ValueError('gamma must be a positive scalar (float or integer)')

    E, V_label, _, _, _ = GetGKInput(G)
    return gkCpy.CalculateHistogramKernelPy(E, V_label, float(gamma), 6)


def CalculateVertexEdgeHistGaussKernel(G, gamma=0.5):
    """Vertex Edge Histogram RBF Kernel"""
    if not isinstance(gamma, (float, int)):
        raise TypeError('gamma must be a positive scalar (float or integer)')

    if gamma <= 0.0:
        raise ValueError('gamma must be a positive scalar (float or integer)')

    E, V_label, _, _, _ = GetGKInput(G)
    return gkCpy.CalculateHistogramKernelPy(E, V_label, float(gamma), 7)


# === Random Walk Kernels ===


def CalculateGeometricRandomWalkKernel(
    G, par=1.0, max_iterations=100, eps=10.0 ** (-10)
):
    """Geometric Random Walk Kernel"""
    if not isinstance(par, (float, int)):
        raise TypeError('par must be a scalar (float or integer)')

    if par == 0:
        warnings.warn('Invoking kernel with par == 0.0')

    if not isinstance(max_iterations, int):
        raise TypeError('max_iterations must be a positive integer')

    if max_iterations <= 0:
        raise ValueError('max_iterations must be a positive integer')

    if not isinstance(eps, (float, int)):
        raise TypeError('eps must be a non-negative scalar (float or integer)')

    if eps < 0.0:
        raise ValueError('eps must be a non-negative scalar (float or integer)')

    E, V_label, _, _, _ = GetGKInput(G)
    return gkCpy.CalculateGeometricRandomWalkKernelPy(
        E, V_label, float(par), max_iterations, float(eps)
    )


def CalculateExponentialRandomWalkKernel(G, par=1.0):
    """Exponential Random Walk Kernel"""
    if not isinstance(par, (float, int)):
        raise TypeError('par must be a scalar (float or integer)')

    if par == 0:
        warnings.warn('Invoking kernel with par == 0.0')

    E, V_label, _, _, _ = GetGKInput(G)
    return gkCpy.CalculateExponentialRandomWalkKernelPy(E, V_label, float(par))


def CalculateKStepRandomWalkKernel(G, par):
    """K-step Random Walk Kernel

    Allow user to provide own list of k-step weights.
    """
    if not isinstance(par, collections.Sequence):
        raise TypeError(
            'par must be a sequence of scalars (floats or integers)'
        )

    gk_par = gkCpy.DoubleVector([float(p) for p in par])
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
