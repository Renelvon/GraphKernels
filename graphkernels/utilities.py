"""
Utility functions for the graph kernels package
"""

import numpy as np

# FIXME: Avoid double-import by exporting names in __init__
from graphkernels import graphkernels as gkCpy


def _get_graph_info(g):
    """Extract graphs information from an igraph object"""

    # matrix of edges
    E = np.zeros(shape=(len(g.es), 2))
    for i in range(len(g.es)):
        E[i, :] = g.es[i].tuple
    # there are multiple edge attributes
    if len(g.es.attributes()) > 1:
        print(
            "There are multiple edge attributes! The first attribute %s is used"
            % g.es.attributes()[0]
        )

    # an edge attribute is missing
    if len(g.es.attributes()) == 0:
        g.es["label"] = 1

    e_attr_name = g.es.attributes()[0]
    e_attr_values = np.asarray(g.es[e_attr_name]).reshape(len(g.es), 1)
    E = np.hstack((E, e_attr_values))

    # if len(g.vs.attributes()) > 1:
    #   print(
    #        "There are multiple vertex attributes! The first attribute %s (or the label) is used"
    #        % g.vs.attributes()[0]
    #    )

    if len(g.vs.attributes()) == 0:
        g.vs["label"] = 1

    # Default to using a 'label' attribute of the graph if it is
    # present. This also accounts for the case where there are 2
    # or more attributes.
    if 'label' in g.vs.attributes():
        v_attr_name = 'label'
    else:
        # FIXME
        # https://github.com/AntoinePrv/GraphKernels/commit/ed097a3680c9e0ee91913dc2d2d4e2efa4a32b32
        v_attr_name = g.vs.attributes()[0]

    v_attr_values = (
        np.asarray(g.vs[v_attr_name]).reshape(len(g.vs), 1).astype(int)
    )

    return E, v_attr_values, len(g.vs), len(g.es), g.maxdegree()


def GetGKInput(G):
    """
    Given a list of graphs, extract graph information.

    Convert in the desired input for graphkernels package.
    """

    E = gkCpy.VecMatrixXi()
    V_label = gkCpy.IntIntVector()
    V_count = gkCpy.IntVector()
    E_count = gkCpy.IntVector()
    D_max = gkCpy.IntVector()

    for graph in G:
        edge, vlabel, vsize, esize, maxdegree = _get_graph_info(graph)
        E.append(edge)
        V_label.append(gkCpy.IntVector(vlabel.reshape(-1).tolist()))
        V_count.append(vsize)
        E_count.append(esize)
        D_max.append(maxdegree)

    return E, V_label, V_count, E_count, D_max


def GetAdjMatList(G):

    adj_mat = gkCpy.VecMatrixXi()
    adj_list = gkCpy.IntIntIntVector()

    for graph in G:
        am_cur = graph.get_adjacency()  # adjacency matrix of i-th graph
        am_cur = np.array(am_cur.data)
        adj_mat.append(am_cur)

        al_cur = np.asarray(graph.get_adjlist())  # adjacency list of i-th graph

        vs = gkCpy.IntIntVector()
        for neighbours in al_cur:
            if isinstance(neighbours, np.ndarray):
                neighbours = neighbours.reshape(-1).tolist()
            vs.append(gkCpy.IntVector(neighbours))
        adj_list.append(vs)

    return adj_mat, adj_list


def normalizekm(K):
    """
    Normalize the kernel matrix.
    Based on a funciton from here
    http://members.cbio.mines-paristech.fr/~nshervashidze/code/
    which was originally Karsten Borgwardt

    Normalize the kernel matrix such that
    diag(result) = 1, i.e. K(x,y) / sqrt(K(x,x) * K(y,y))

    K:
        a kernel matrix

    returns:
        the normalized result
    """
    nv = np.sqrt(np.diag(K))
    nm = nv[:, np.newaxis] * nv[:, np.newaxis].T
    Knm = nm ** -1

    Knm[np.where(np.isnan(Knm))] = 0

    return K * Knm
