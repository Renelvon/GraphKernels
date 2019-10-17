"""
Demo: Calculate the kernel matrices on the MUTAG data.
"""

from os import path

from igraph import Graph

from graphkernels import kernels as gk

N_GRAPHS = 188

MUTAG_DIR = 'mutag'
THIS_DIR = path.dirname(__file__)


def load_mutag_graphs():
    paths = (
        path.join(MUTAG_DIR, ('mutag_%d.graphml' % i))
        for i in range(1, N_GRAPHS + 1)
    )
    return tuple(Graph.Read_GraphML(path) for path in paths)


def compute_all_kernels(graphs):
    return (
        gk.CalculateEdgeHistKernel(graphs),
        gk.CalculateVertexHistKernel(graphs),
        gk.CalculateVertexEdgeHistKernel(graphs),
        gk.CalculateVertexVertexEdgeHistKernel(graphs),
        gk.CalculateEdgeHistGaussKernel(graphs),
        gk.CalculateVertexHistGaussKernel(graphs),
        gk.CalculateVertexEdgeHistGaussKernel(graphs),
#        gk.CalculateGeometricRandomWalkKernel(graphs),
#        gk.CalculateExponentialRandomWalkKernel(graphs),
#        gk.CalculateKStepRandomWalkKernel(graphs),
        gk.CalculateWLKernel(graphs),
        gk.CalculateConnectedGraphletKernel(graphs, 4),
        gk.CalculateGraphletKernel(graphs, 4),
#        gk.CalculateShortestPathKernel(graphs)
    )


def main():
    """
    Compute all kernels.
    """
    graphs = load_mutag_graphs()
    for result in  compute_all_kernels(graphs):
        print(result.sum())


if __name__ == '__main__':
    main()
