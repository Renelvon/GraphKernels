/* Copyright (c) 2018-2019 the `graphkernels` developers
 * All rigths reserved.
 */

#include "graphkernels.h"

// ===== connected graphlet kernel for k = 5 ===== //
void getMinValue(MatrixXi& iam, vector<int>& idx, vector<int>& sums) {
  SparseMatrix<int> am;
  am = iam.sparseView();

  sums.clear();
  sums.resize(idx.size());
  fill(sums.begin(), sums.end(), 0);
  for (auto i = 0; i < idx.size(); ++i) {
    Int k = idx[i];
    for (SparseMatrix<int>::InnerIterator it(am, k); it; ++it) {
      if (find(idx.begin(), idx.end(), it.row()) != idx.end()) {
        sums[i] += it.value();
      }
    }
  }
  sums.push_back(1);
}

VectorXd countConnectedGraphletsFive(MatrixXi& am,
                                     vector<vector<int>>& al,
                                     VectorXd& count_gr) {
  auto n = int{al.size()};

  vector<double> w = {
      1.0 / 120.0, 1.0 / 72.0, 1.0 / 48.0, 1.0 / 36.0, 1.0 / 28.0, 1.0 / 20.0,
      1.0 / 14.0,  1.0 / 10.0, 1.0 / 12.0, 1.0 / 8.0,  1.0 / 8.0,  1.0 / 4.0,
      1.0 / 2.0,   1.0 / 12.0, 1.0 / 12.0, 1.0 / 4.0,  1.0 / 4.0,  1.0 / 2.0,
      0.0,         0.0,        0.0};
  // int n = (int)am.rows();
  vector<int> L1(n);
  iota(L1.begin(), L1.end(), 0);
  vector<int> idx(5);
  vector<int> sums;

  for (auto&& i : L1) {
    for (auto&& j : al[i]) {
      for (auto&& k : al[j]) {
        if (k != i) {
          for (auto&& l : al[k]) {
            if (l != i && l != j) {
              for (auto&& m : al[l]) {
                if (m != i && m != j && m != k) {
                  int aux = am.coeff(i, k) + am.coeff(i, l) + am.coeff(i, m) +
                            am.coeff(j, l) + am.coeff(j, m) + am.coeff(k, m);
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
                    int aux1 = *min_element(sums.begin(), sums.end());
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
                    copy(sums.begin(), sums.end(), back_inserter(aux1));
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
                          if (sums[ii] == 3)
                            ind.push_back(ii);
                        }
                        // idx[0] = i; idx[1] = j; idx[2] = k; idx[3] = l;
                        // idx[4] = m;
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
                    copy(sums.begin(), sums.end(), back_inserter(aux1));
                    sort(aux1.begin(), aux1.end());
                    if (aux1[0] == 2) {
                      count_gr[7] += w[7];
                    } else if (aux1[1] == 1) {
                      count_gr[17] += w[17];
                    } else {
                      vector<int> ind;
                      for (auto ii = 0; ii < sums.size(); ++ii) {
                        if (sums[ii] == 3)
                          ind.push_back(ii);
                      }
                      for (auto ii = 0; ii < sums.size(); ++ii) {
                        if (sums[ii] == 1)
                          ind.push_back(ii);
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
                if (m != i && m != j && m != l && am.coeff(i, m) == 0 &&
                    am.coeff(j, m) == 0 && am.coeff(l, m) == 0) {
                  count_gr[19] += w[19];
                }
              }
            }
          }
        }
      }
    }
    // count graphlets of type 19 and 21
  for (auto m = 3 - 0; m < al[i].size(); ++m) {
    for (auto l = m - 1; l < al[i].size(); ++l) {
      for (auto k = l - 1; k < al[i].size(); ++k) {
        for (auto j = k - 1; j < al[i].size(); ++j) {
            auto aux =
                am.coeff(al[i][j], al[i][k]) + am.coeff(al[i][j], al[i][l]) +
                am.coeff(al[i][j], al[i][m]) + am.coeff(al[i][k], al[i][l]) +
                am.coeff(al[i][k], al[i][m]) + am.coeff(al[i][l], al[i][m]);
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
  return count_gr;
}

// ===== connected graphlet kernel for k = 4 ===== //
VectorXd countConnectedGraphletsFour(MatrixXi& am,
                                     vector<vector<int>>& al,
                                     VectorXd& count_gr) {
  vector<double> w = {1.0 / 24.0, 1.0 / 12.0, 1.0 / 4.0,
                      0.0,        1.0 / 8.0,  1.0 / 2.0};
  auto n = int{am.rows()};
  vector<int> L1(n);
  iota(L1.begin(), L1.end(), 0);

  for (auto&& i : L1) {
    for (auto&& j : al[i]) {
      for (auto&& k : al[j]) {
        if (k != i) {
          for (auto&& l : al[k]) {
            if (l != i && l != j) {
              int aux = am.coeff(i, k) + am.coeff(i, l) + am.coeff(j, l);
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
    for (auto l = 2 - 0; l < al[i].size(); ++l) {
      for (auto k = l - 1; k < al[i].size(); ++k) {
        for (auto j = k - 1; j < al[i].size(); ++j) {
          if (am.coeff(al[i][j], al[i][k]) == 0 &&
              am.coeff(al[i][j], al[i][l]) == 0 &&
              am.coeff(al[i][k], al[i][l]) == 0) {
            count_gr[3]++;
          }
        }
      }
    }
  }
  return count_gr;
}

// ===== connected graphlet kernel for k = 3 ===== //
VectorXd countConnectedGraphletsThree(MatrixXi& am,
                                      vector<vector<int>>& al,
                                      VectorXd& count_gr) {
  vector<double> w = {1.0 / 2.0, 1.0 / 6.0};
  auto n = int{am.rows()};
  vector<int> L1(n);
  iota(L1.begin(), L1.end(), 0);

  for (auto&& i : L1) {
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
  return count_gr;
}

// Python function
MatrixXd CalculateConnectedGraphletKernelPy(
    vector<MatrixXi>& graph_adj_all,
    vector<vector<vector<int>>>& graph_adjlist_all,
    int k) {
  // decrement one to start indices from zero
  // for (auto&& X : graph_adjlist_all)
  // for (auto&& vec : X)
  // for (auto&& x : vec) x--;

  // MatrixXd freq;
  int freq_size = 0;

  switch (k) {
    case 3:
      freq_size = 2;
      break;
    case 4:
      freq_size = 6;
      break;
    case 5:
      freq_size = 21;
      break;
  }

  MatrixXd freq(graph_adjlist_all.size(), freq_size);
  // freq = MatrixXd::Zero(graph_adjlist_all.size(), freq_size);

  vector<int> idx_graph(graph_adjlist_all.size());
  iota(idx_graph.begin(), idx_graph.end(), 0);

  VectorXd count_g;
  VectorXd freq_row;

  for (auto&& i : idx_graph) {
    // VectorXd count_g;
    // VectorXd freq_row;
    // Eigen::Map<VectorXd>(freq_row.data(), freq.row(i).size()) = freq.row(i);
    freq_row = VectorXd::Zero(freq_size);

    if (k == 3) {
      count_g = countConnectedGraphletsThree(graph_adj_all[i],
                                             graph_adjlist_all[i], freq_row);

    } else if (k == 4) {
      count_g = countConnectedGraphletsFour(graph_adj_all[i],
                                            graph_adjlist_all[i], freq_row);

    } else if (k == 5) {
      count_g = countConnectedGraphletsFive(graph_adj_all[i],
                                            graph_adjlist_all[i], freq_row);
    }

    // freq.row(i) = Eigen::Map<VectorXd>(count_g.data(), count_g.size());
    freq.row(i) = count_g;

    if (freq.row(i).sum() != 0) {
      freq.row(i) /= freq.row(i).sum();
    }
  }
  MatrixXd K = freq * freq.transpose();

  return K;
}
