/* Copyright (c) 2018-2019 the `graphkernels` developers
 * All rights reserved.
 */

#include "connected_graphlet.h"

#include <Eigen/Sparse>

#include <algorithm>

using std::sort;
using std::vector;

using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::SparseMatrix;
using Eigen::VectorXd;

constexpr auto FREQ_SIZE_3 = 2;
constexpr auto FREQ_SIZE_4 = 6;
constexpr auto FREQ_SIZE_5 = 21;

void getMinValue(const MatrixXi& iam, const vector<int>& idx, vector<int>& sums) {
    SparseMatrix<int> am = iam.sparseView();

    sums.clear();
    sums.resize(idx.size());
    fill(sums.begin(), sums.end(), 0);
    for (auto i = 0; i < idx.size(); ++i) {
        for (SparseMatrix<int>::InnerIterator it(am, idx[i]); it; ++it) {
            if (find(idx.cbegin(), idx.cend(), it.row()) != idx.cend()) {
                sums[i] += it.value();
            }
        }
    }
    sums.push_back(1);
}

VectorXd countConnectedGraphletsFive(
        const MatrixXi& am,
        const vector<vector<int>>& al) {
    vector<double> w = {
        1.0 / 120.0, 1.0 / 72.0, 1.0 / 48.0, 1.0 / 36.0, 1.0 / 28.0, 1.0 / 20.0,
        1.0 / 14.0,  1.0 / 10.0, 1.0 / 12.0, 1.0 / 8.0,  1.0 / 8.0,  1.0 / 4.0,
        1.0 / 2.0,   1.0 / 12.0, 1.0 / 12.0, 1.0 / 4.0,  1.0 / 4.0,  1.0 / 2.0,
        0.0,         0.0,        0.0};

    const auto n = al.size();
    VectorXd count_gr = VectorXd::Zero(FREQ_SIZE_5);

    vector<int> idx(5);
    vector<int> sums;

    for (auto i = 0; i < n; ++i) {
        for (auto&& j : al[i]) {
            for (auto&& k : al[j]) {
                if (k != i) {
                    for (auto&& l : al[k]) {
                        if (l != i && l != j) {
                            for (auto&& m : al[l]) {
                                if (m != i && m != j && m != k) {
                                    const auto aux =
                                        am.coeff(i, k) + am.coeff(i, l) +
                                        am.coeff(i, m) + am.coeff(j, l) +
                                        am.coeff(j, m) + am.coeff(k, m);
                                    if (aux == 6) {
                                        count_gr[0] += w[0];
                                    } else if (aux == 5) {
                                        count_gr[1] += w[1];
                                    } else if (aux == 4) {
                                        idx[0] = i;
                                        idx[1] = j;
                                        idx[2] = k;
                                        idx[3] = l;
                                        idx[4] = m;
                                        getMinValue(am, idx, sums);
                                        const auto aux1 = *min_element(
                                                sums.cbegin(), sums.cend());
                                        if (aux1 == 2) {
                                            count_gr[3] += w[3];
                                        } else {
                                            count_gr[2] += w[2];
                                        }
                                    } else if (aux == 3) {
                                        idx[0] = i;
                                        idx[1] = j;
                                        idx[2] = k;
                                        idx[3] = l;
                                        idx[4] = m;
                                        getMinValue(am, idx, sums);
                                        sort(sums.begin(), sums.end());
                                        if (sums[0] == 1) {
                                            count_gr[8] += w[8];
                                        } else if (sums[1] == 3) {
                                            count_gr[4] += w[4];
                                        } else if (sums[2] == 2) {
                                            count_gr[13] += w[13];
                                        } else {
                                            count_gr[5] += w[5];
                                        }
                                    } else if (aux == 2) {
                                        idx[0] = i;
                                        idx[1] = j;
                                        idx[2] = k;
                                        idx[3] = l;
                                        idx[4] = m;
                                        getMinValue(am, idx, sums);
                                        vector<int> aux1;
                                        copy(sums.cbegin(), sums.cend(),
                                                back_inserter(aux1));
                                        sort(aux1.begin(), aux1.end());
                                        if (aux1[0] == 1) {
                                            if (aux1[2] == 2) {
                                                count_gr[15] += w[15];
                                            } else {
                                                count_gr[9] += w[9];
                                            }
                                        } else {
                                            if (aux1[3] == 2) {
                                                count_gr[10] += w[10];
                                            } else {
                                                vector<int> ind;
                                                for (auto ii = 0; ii < sums.size(); ++ii) {
                                                    if (sums[ii] == 3) {
                                                        ind.push_back(ii);
                                                    }
                                                }
                                                if (am.coeff(idx[ind[0]], idx[ind[1]]) == 1) {
                                                    count_gr[6] += w[6];
                                                } else {
                                                    count_gr[14] += w[14];
                                                }
                                            }
                                        }
                                    } else if (aux == 1) {
                                        idx[0] = i;
                                        idx[1] = j;
                                        idx[2] = k;
                                        idx[3] = l;
                                        idx[4] = m;
                                        getMinValue(am, idx, sums);
                                        vector<int> aux1;
                                        copy(sums.cbegin(), sums.cend(),
                                                back_inserter(aux1));
                                        sort(aux1.begin(), aux1.end());
                                        if (aux1[0] == 2) {
                                            count_gr[7] += w[7];
                                        } else if (aux1[1] == 1) {
                                            count_gr[17] += w[17];
                                        } else {
                                            vector<int> ind;
                                            for (auto ii = 0; ii < sums.size(); ++ii) {
                                                if (sums[ii] == 3) {
                                                    ind.push_back(ii);
                                                }
                                            }
                                            for (auto ii = 0; ii < sums.size(); ++ii) {
                                                if (sums[ii] == 1) {
                                                    ind.push_back(ii);
                                                }
                                            }
                                            if (am.coeff(idx[ind[0]], idx[ind[1]]) == 1) {
                                                count_gr[16] += w[16];
                                            } else {
                                                count_gr[11] += w[11];
                                            }
                                        }
                                    } else {
                                        count_gr[12] += w[12];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        // count graphlets of type 20
        for (auto&& j : al[i]) {
            for (auto&& k : al[j]) {
                if (k != i && am.coeff(i, k) == 0) {
                    for (auto&& l : al[k]) {
                        if (l != i && l != j && am.coeff(i, l) == 0 &&
                                am.coeff(j, l) == 0) {
                            for (auto&& m : al[k]) {
                                if (m != i && m != j && m != l
                                        && am.coeff(i, m) == 0
                                        && am.coeff(j, m) == 0
                                        && am.coeff(l, m) == 0) {
                                    count_gr[19] += w[19];
                                }
                            }
                        }
                    }
                }
            }
        }
        // count graphlets of type 19 and 21
        for (auto m = al[i].size(); m-- > 3;) {
            for (auto l = m; l-- > 2;) {
                for (auto k = l; k-- > 1;) {
                    for (auto j = k; j-- > 0;) {
                        auto aux =
                            am.coeff(al[i][j], al[i][k]) +
                            am.coeff(al[i][j], al[i][l]) +
                            am.coeff(al[i][j], al[i][m]) +
                            am.coeff(al[i][k], al[i][l]) +
                            am.coeff(al[i][k], al[i][m]) +
                            am.coeff(al[i][l], al[i][m]);
                        if (aux == 1) {
                            count_gr[18]++;
                        } else if (aux == 0) {
                            count_gr[20]++;
                        }
                    }
                }
            }
        }
    }
    const auto csum = count_gr.sum();
    return (csum == 0.0) ? count_gr : count_gr / csum;
}

VectorXd countConnectedGraphletsFour(
        const MatrixXi& am,
        const vector<vector<int>>& al) {
    vector<double> w = {1.0 / 24.0, 1.0 / 12.0, 1.0 / 4.0,
                        0.0,        1.0 / 8.0,  1.0 / 2.0};

    VectorXd count_gr = VectorXd::Zero(FREQ_SIZE_4);
    const auto n = am.rows();
    for (auto i = 0; i < n; ++i) {
        for (auto&& j : al[i]) {
            for (auto&& k : al[j]) {
                if (k != i) {
                    for (auto&& l : al[k]) {
                        if (l != i && l != j) {
                            const auto aux =
                                am.coeff(i, k) +
                                am.coeff(i, l) +
                                am.coeff(j, l);
                            if (aux == 3) {
                                count_gr[0] += w[0];
                            } else if (aux == 2) {
                                count_gr[1] += w[1];
                            } else if (aux == 1) {
                                if (am.coeff(i, l) == 1) {
                                    count_gr[4] += w[4];
                                } else {
                                    count_gr[2] += w[2];
                                }
                            } else {
                                count_gr[5] += w[5];
                            }
                        }
                    }
                }
            }
        }

        // count "stars"
        for (auto l = al[i].size(); l-- > 2;) {
            for (auto k = l; k-- > 1;) {
                for (auto j = k; j-- > 0;) {
                    if (am.coeff(al[i][j], al[i][k]) == 0 &&
                            am.coeff(al[i][j], al[i][l]) == 0 &&
                            am.coeff(al[i][k], al[i][l]) == 0) {
                        count_gr[3]++;
                    }
                }
            }
        }
    }
    const auto csum = count_gr.sum();
    return (csum == 0.0) ? count_gr : count_gr / csum;
}

VectorXd countConnectedGraphletsThree(
        const MatrixXi& am,
        const vector<vector<int>>& al) {
    vector<double> w = {1.0 / 2.0, 1.0 / 6.0};

    VectorXd count_gr = VectorXd::Zero(FREQ_SIZE_3);
    const auto n = am.rows();
    for (auto i = 0; i < n; ++i) {
        for (auto&& j : al[i]) {
            for (auto&& k : al[j]) {
                if (k != i) {
                    if (am.coeff(i, k) == 1) {
                        count_gr[1] += w[1];
                    } else {
                        count_gr[0] += w[0];
                    }
                }
            }
        }
    }
    const auto csum = count_gr.sum();
    return (csum == 0.0) ? count_gr : count_gr / csum;
}

MatrixXd CalculateConnectedGraphletKernelThreePy(
        const vector<MatrixXi>& graph_adj_all,
        const vector<vector<vector<int>>>& graph_adjlist_all) {
    MatrixXd freq(FREQ_SIZE_3, graph_adjlist_all.size());

    for (auto i = 0; i < graph_adjlist_all.size(); ++i) {
        freq.col(i) = countConnectedGraphletsThree(
                graph_adj_all[i], graph_adjlist_all[i]);
    }

    return freq.transpose() * freq;
}

MatrixXd CalculateConnectedGraphletKernelFourPy(
        const vector<MatrixXi>& graph_adj_all,
        const vector<vector<vector<int>>>& graph_adjlist_all) {
    MatrixXd freq(FREQ_SIZE_4, graph_adjlist_all.size());

    for (auto i = 0; i < graph_adjlist_all.size(); ++i) {
        freq.col(i) = countConnectedGraphletsFour(
                graph_adj_all[i], graph_adjlist_all[i]);
    }

    return freq.transpose() * freq;
}

MatrixXd CalculateConnectedGraphletKernelFivePy(
        const vector<MatrixXi>& graph_adj_all,
        const vector<vector<vector<int>>>& graph_adjlist_all) {
    MatrixXd freq(FREQ_SIZE_5, graph_adjlist_all.size());

    for (auto i = 0; i < graph_adjlist_all.size(); ++i) {
        freq.col(i) = countConnectedGraphletsFive(
                graph_adj_all[i], graph_adjlist_all[i]);
    }

    return freq.transpose() * freq;
}
