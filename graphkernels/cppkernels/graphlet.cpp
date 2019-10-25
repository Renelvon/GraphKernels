/* Copyright (c) 2018-2019 the `graphkernels` developers
 * All rigths reserved.
 */

#include "graphkernels.h"


// ========================================================= //
// ==================== Graphlet kernel ==================== //
// ========================================================= //
// ===== graphlet kernel for k = 4 ===== //
int find_min(int a, int b, int c) {
  int m;
  int mini = a;
  if (b < mini)
    mini = b;
  if (c < mini)
    mini = c;
  if (mini == a) {
    if (mini == b) {
      if (mini == c) {
        m = 7;
      } else {
        m = 4;
      }
    } else {
      if (mini == c) {
        m = 5;
      } else {
        m = 1;
      }
    }
  } else {
    if (mini == b) {
      if (mini == c) {
        m = 6;
      } else {
        m = 2;
      }
    } else {
      m = 3;
    }
  }
  return m;
}

void card_ThreeInter(vector<int>& L1,
                     vector<int>& L2,
                     vector<int>& L3,
                     vector<int>& card) {
  card.resize(7);
  fill(card.begin(), card.end(), 0);
  int i = 0, j = 0, k = 0;

  while (i < (int)L1.size() && j < (int)L2.size() && k < (int)L3.size()) {
    int m = find_min(L1[i], L2[j], L3[k]);
    card[m - 1] += 1;
    switch (m) {
      case 1:
        i++;
        break;
      case 2:
        j++;
        break;
      case 3:
        k++;
        break;
      case 4:
        i++;
        j++;
        break;
      case 5:
        i++;
        k++;
        break;
      case 6:
        j++;
        k++;
        break;
      case 7:
        i++;
        j++;
        k++;
        break;
    }
  }

  if (i < (int)L1.size() || j < (int)L2.size() || k < (int)L3.size()) {
    if (i >= (int)L1.size() && j >= (int)L2.size()) {
      card[2] += (int)L3.size() - k;
      k = (int)L3.size();
    } else {
      if (i >= (int)L1.size() && k >= (int)L3.size()) {
        card[1] += (int)L2.size() - j;
        j = (int)L2.size();
      } else {
        if (j >= (int)L2.size() && k >= (int)L3.size()) {
          card[0] += (int)L1.size() - i;
          i = (int)L1.size();
        } else {
          if (i >= (int)L1.size()) {
            while (j < (int)L2.size() && k < (int)L3.size()) {
              if (L2[j] < L3[k]) {
                card[1]++;
                j++;
              } else {
                if (L2[j] > L3[k]) {
                  card[2]++;
                  k++;
                } else {
                  card[5]++;
                  j++;
                  k++;
                }
              }
            }
          } else {
            if (j >= (int)L2.size()) {
              while (i < (int)L1.size() && k < (int)L3.size()) {
                if (L1[i] < L3[k]) {
                  card[0]++;
                  i++;
                } else {
                  if (L1[i] > L3[k]) {
                    card[2]++;
                    k++;
                  } else {
                    card[4]++;
                    i++;
                    k++;
                  }
                }
              }
            } else {
              if (k >= (int)L3.size()) {
                while (i < (int)L1.size() && j < (int)L2.size()) {
                  if (L1[i] < L2[j]) {
                    card[0]++;
                    i++;
                  } else {
                    if (L1[i] > L2[j]) {
                      card[1]++;
                      j++;
                    } else {
                      card[3]++;
                      i++;
                      j++;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (i < (int)L1.size() || j < (int)L2.size() || k < (int)L3.size()) {
    if (i >= (int)L1.size() && j >= (int)L2.size()) {
      card[2] += (int)L3.size() - k;
    } else if (i >= (int)L1.size() && k >= (int)L3.size()) {
      card[1] += (int)L2.size() - j;
    } else if (j >= (int)L2.size() && k >= (int)L3.size()) {
      card[0] += (int)L1.size() - i;
    }
  }
}

void getIndices(vector<int>& o_set1,
                vector<int>& o_set2,
                vector<int>& inter,
                vector<int>& diff1,
                vector<int>& diff2) {
  vector<int> inter_(min(o_set1.size(), o_set2.size()), -1);
  vector<int> diff1_(max(o_set1.size(), o_set2.size()), -1);
  vector<int> diff2_(max(o_set1.size(), o_set2.size()), -1);

  int i = 0, j = 0;
  while (i < (int)o_set1.size() && j < (int)o_set2.size()) {
    if (o_set1[i] < o_set2[j]) {
      diff1_[i] = o_set1[i];
      i++;
    } else if (o_set1[i] > o_set2[j]) {
      diff2_[j] = o_set2[j];
      j++;
    } else {
      inter_[i] = o_set1[i];
      i++;
      j++;
    }
  }

  if (i < (int)o_set1.size()) {
    for (int k = i; k < (int)o_set1.size(); ++k) {
      diff1_[k] = o_set1[k];
    }
  } else if (j < (int)o_set2.size()) {
    for (int k = j; k < (int)o_set2.size(); ++k) {
      diff2_[k] = o_set2[k];
    }
  }

  inter.clear();
  for (auto&& x : inter_) {
    if (x >= 0)
      inter.push_back(x);
  }
  diff1.clear();
  for (auto&& x : diff1_) {
    if (x >= 0)
      diff1.push_back(x);
  }
  diff2.clear();
  for (auto&& x : diff2_) {
    if (x >= 0)
      diff2.push_back(x);
  }
}

VectorXd countGraphletsFour(vector<vector<int>>& al, VectorXd& count_gr) {
  double n = (double)al.size();
  vector<double> w = {1.0 / 12.0, 1.0 / 10.0, 1.0 / 8.0, 1.0 / 6.0,
                      1.0 / 8.0,  1.0 / 6.0,  1.0 / 6.0, 1.0 / 4.0,
                      1.0 / 4.0,  1.0 / 2.0,  0};
  vector<int> inter, diff1, diff2, card;
  vector<double> inter_count(11);
  vector<int> v;
  vector<int>::iterator it;

  double m = 0.0;
  for (auto&& vec : al) {
    m += (double)vec.size();
  }
  m /= 2.0;

  vector<int> v1(al.size());
  iota(v1.begin(), v1.end(), 0);
  for (auto&& i : v1) {
    for (auto&& j : al[i]) {
      double K = 0.0;
      fill(inter_count.begin(), inter_count.end(), 0.0);
      getIndices(al[i], al[j], inter, diff1, diff2);
      for (auto&& k : inter) {
        card_ThreeInter(al[i], al[j], al[k], card);
        inter_count[0] += 0.5 * (double)card[6];
        inter_count[1] += 0.5 * (double)(card[3] - 1.0);
        inter_count[1] += 0.5 * (double)(card[4] - 1.0);
        inter_count[1] += 0.5 * (double)(card[5] - 1.0);
        inter_count[2] += 0.5 * (double)card[0];
        inter_count[2] += 0.5 * (double)card[1];
        inter_count[2] += (double)card[2];
        inter_count[6] += n - (double)accumulate(card.begin(), card.end(), 0);
        K += 0.5 * (double)card[6] + 0.5 * (double)(card[4] - 1.0) +
             0.5 * (double)(card[5] - 1.0) + card[2];
      }
      v.clear();
      v.resize(diff1.size());
      sort(diff1.begin(), diff1.end());
      sort(al[i].begin(), al[i].end());
      it = set_difference(diff1.begin(), diff1.end(), al[i].begin(),
                          al[i].end(), v.begin());
      v.resize(it - v.begin());
      for (auto&& k : v) {
        card_ThreeInter(al[i], al[j], al[k], card);
        inter_count[1] += 0.5 * (double)card[6];
        inter_count[2] += 0.5 * (double)card[3];
        inter_count[2] += 0.5 * (double)card[4];
        inter_count[4] += 0.5 * (double)(card[5] - 1.0);
        inter_count[3] += 0.5 * (double)(card[0] - 2.0);
        inter_count[5] += 0.5 * (double)card[1];
        inter_count[5] += (double)card[2];
        inter_count[7] += n - (double)accumulate(card.begin(), card.end(), 0);
        K += 0.5 * (double)card[6] + 0.5 * (double)card[4] +
             0.5 * (double)(card[5] - 1.0) + card[2];
      }
      v.clear();
      v.resize(diff2.size());
      sort(diff2.begin(), diff2.end());
      it = set_difference(diff2.begin(), diff2.end(), v1.begin(), v1.end(),
                          v.begin());
      v.resize(it - v.begin());
      for (auto&& k : v) {
        card_ThreeInter(al[i], al[j], al[k], card);
        inter_count[1] += 0.5 * (double)card[6];
        inter_count[2] += 0.5 * (double)card[3];
        inter_count[4] += 0.5 * (double)(card[4] - 1.0);
        inter_count[2] += 0.5 * (double)card[5];
        inter_count[5] += 0.5 * (double)card[0];
        inter_count[3] += 0.5 * (double)(card[1] - 2.0);
        inter_count[5] += (double)card[2];
        inter_count[7] += n - (double)accumulate(card.begin(), card.end(), 0);
        K += 0.5 * (double)card[6] + 0.5 * (double)(card[4] - 1.0) +
             0.5 * (double)card[5] + card[2];
      }
      inter_count[8] += m + 1.0 - (double)v1.size() - (double)al[i].size() - K;
      inter_count[9] +=
          (n - (double)inter.size() - (double)diff1.size() -
           (double)diff2.size()) *
              (n - (double)inter.size() - (double)diff1.size() -
               (double)diff2.size() - 1.0) /
              2 -
          (m + 1.0 - (double)v1.size() - (double)al[i].size() - K);

      for (int k = 0; k < (int)count_gr.size(); ++k) {
        count_gr(k) += inter_count[k] * w[k];
      }
    }
  }

  count_gr(10) = n * (n - 1.0) * (n - 2.0) * (n - 3.0) / (4.0 * 3.0 * 2.0) -
                 count_gr.head(10).sum();
  return count_gr;
}

// ===== graphlet kernel for k = 3 ===== //
void getCardinality(vector<int>& o_set1,
                    vector<int>& o_set2,
                    vector<double>& card) {
  card.resize(3);
  fill(card.begin(), card.end(), 0.0);
  int i = 0, j = 0;
  while (i < (int)o_set1.size() && j < (int)o_set2.size()) {
    if (o_set1[i] < o_set2[j]) {
      card[0] += 1.0;
      i++;
    } else if (o_set1[i] > o_set2[j]) {
      card[1] += 1.0;
      j++;
    } else {
      i++;
      j++;
      card[2] += 1.0;
    }
  }
  card[0] += (double)((int)o_set1.size() - i);
  card[1] += (double)((int)o_set2.size() - j);
}

VectorXd countGraphletsThree(vector<vector<int>>& al, VectorXd& count_gr) {
  double n = (double)al.size();
  vector<double> w = {1.0 / 6.0, 1.0 / 4.0, 1.0 / 2.0};
  vector<double> card(3);

  vector<int> L1(al.size());
  iota(L1.begin(), L1.end(), 0);
  for (auto&& i : L1) {
    for (auto&& j : al[i]) {
      getCardinality(al[i], al[j], card);
      count_gr(0) += w[0] * card[2];
      count_gr(1) += w[1] * (card[0] + card[1] - 2.0);
      count_gr(2) += w[2] * (n - accumulate(card.begin(), card.end(), 0.0));
    }
  }
  count_gr(3) = n * (n - 1.0) * (n - 2.0) / 6.0 -
                (count_gr(0) + count_gr(1) + count_gr(2));
  return count_gr;
}

// Python function
MatrixXd CalculateGraphletKernelPy(
    vector<MatrixXi>& graph_adj_all,
    vector<vector<vector<int>>>& graph_adjlist_all,
    int k) {
  // decrement one to start indices from zero
  // for (auto&& X : graph_adjlist_all)
  // for (auto&& vec : X)
  // for (auto&& x : vec) x--;

  int freq_size = 0;

  switch (k) {
    case 3:
      freq_size = 4;
      break;
    case 4:
      freq_size = 11;
      break;
  }

  // freq = MatrixXd::Zero(graph_adjlist_all.size(), freq_size);
  MatrixXd freq(graph_adjlist_all.size(), freq_size);

  vector<int> idx_graph(graph_adjlist_all.size());
  iota(idx_graph.begin(), idx_graph.end(), 0);

  VectorXd count_g;
  VectorXd freq_row;

  for (auto&& i : idx_graph) {
    freq_row = VectorXd::Zero(freq_size);

    // Eigen::Map<VectorXd>(freq_row.data(), freq.row(i).size()) = freq.row(i);

    if (k == 3) {
      count_g = countGraphletsThree(graph_adjlist_all[i], freq_row);
    } else if (k == 4) {
      count_g = countGraphletsFour(graph_adjlist_all[i], freq_row);
    }

    freq.row(i) = count_g;

    if (freq.row(i).sum() != 0) {
      freq.row(i) /= freq.row(i).sum();
    }
  }

  MatrixXd K = freq * freq.transpose();

  return K;
}
