"""
Functions for computing the graph kernels
"""

import numpy as np
from igraph import Graph

# FIXME: Avoid double-import by exporting names in __init__
from GKextCPy import GKextCPy as gkCpy

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
    """Weisfeiler-Lehman Kernel

    Parametres
    ----------
    par : number of WL iterations
    """
    if not isinstance(par, int):
        raise TypeError('Number of WL iterations must be an integer')

    if par < 0:
        raise ValueError('Number of WL iterations must be non-negative')

    E, V_label, V_count, E_count, D_max = GetGKInput(G) # Extract graph info.
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

    adj_mat, adj_list = GetAdjMatList(G) # Extract graph info.
    return gkCpy.CalculateGraphletKernelPy(adj_mat, adj_list, par)


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

    adj_mat, adj_list = GetAdjMatList(G) # Extract graph info.
    return gkCpy.CalculateConnectedGraphletKernelPy(adj_mat, adj_list, par)


def _floyd_transform(g):
    am = g.shortest_paths_dijkstra()
    am = np.asarray(am).reshape(len(am), len(am))
    am_pos_l = (am > 0).to_list()

    g_floyd = Graph.Adjacency(am_pos_l)

    g_floyd.es['label'] = am[am.nonzero()]
    g_floyd.vs['id'] = np.arange(len(g.vs['label']))
    g_floyd.vs['label'] = g.vs['label']

    return g_floyd


def CalculateShortestPathKernel(G):
    """Shortest Path Kernel"""

    floyd_graphs = tuple(_floyd_transform(g) for g in G)
    G_floyd = np.array(floyd_graphs)

    return CalculateKStepRandomWalkKernel(G_floyd, par=(0, 1))
