"""
Converting mutag to list format
"""

import numpy as np
from igraph import Graph

N_GRAPHS = 188

def main():
    mutag_list = [
        Graph.Read_GraphML(
            "/home/eghisu/projects/gk_python_wrapper/gk_python_c/data/mutag/mutag_%d.graphml" % i
        ) for i in range(1, N_GRAPHS + 1)
    ]
    np.save(
        "/home/eghisu/projects/gk_python_wrapper/gk_python_c/data/mutag_pydata",
        mutag_list
    )


if __name__ == '__main__':
    main()
