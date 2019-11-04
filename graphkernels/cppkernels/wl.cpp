/* Copyright (c) 2018-2019 the `graphkernels` developers
 * All rights reserved.
 */

#include "wl.h"

#include <algorithm>
#include <functional>
#include <numeric>
#include <set>

using std::accumulate;
using std::greater;
using std::set;
using std::sort;
using std::vector;

using Eigen::MatrixXd;
using Eigen::MatrixXi;

// bucket sort used in Weisfeiler-Leiman graph kernel
void bucketsort(vector<int>& x, vector<int>& index, int label_max) {
    vector<vector<int>> buckets;
    buckets.resize(label_max + 1);

    for (auto itr = index.begin(), end = index.end(); itr != end; ++itr) {
        buckets[x[*itr]].push_back(*itr);
    }

    int counter = 0;
    for (auto itr = buckets.begin(), end = buckets.end(); itr != end; ++itr) {
        for (auto itr2 = (*itr).begin(), end2 = (*itr).end(); itr2 != end2;
                ++itr2) {
            index[counter] = *itr2;
            counter++;
        }
    }
}

// Weisfeiler-Leiman graph kernel
MatrixXd WLKernelMatrix(
        vector<MatrixXi>& E,
        vector<vector<int>>& V_label,
        vector<int>& num_v,
        vector<int>& num_e,
        vector<int>& degree_max,
        int h_max) {
    MatrixXd K_mat(V_label.size(), V_label.size());

    K_mat.setZero();

    auto n = static_cast<int>(E.size());
    int v_all = accumulate(num_v.begin(), num_v.end(), 0);
    int degree_max_all = *max_element(degree_max.begin(), degree_max.end());
    vector<int> label_max_vec(n);
    for (int i = 0; i < n; i++) {
        label_max_vec[i] = *max_element(V_label[i].begin(), V_label[i].end());
    }
    int label_max = *max_element(label_max_vec.begin(), label_max_vec.end());

    int raise = 0;
    vector<int> counter(*max_element(num_v.begin(), num_v.end()));
    MatrixXi nei_list(v_all, degree_max_all + 1);
    MatrixXi label_list(v_all, h_max + 1);
    vector<int> x(v_all);
    vector<int> index(v_all);
    vector<int> index_org(v_all);
    vector<int> graph_index(v_all);

    // NEWadd
    label_list.setZero();

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < num_v[i]; j++) {
            label_list(j + raise, 0) = V_label[i][j];
            graph_index[j + raise] = i;
        }
        raise += num_v[i];
    }

    // ===== Increment kernel values using the initial vertex labels =====
    // radix sort
    for (int i = 0; i < v_all; i++) {
        index[i] = i;
        index_org[i] = i;
        x[i] = label_list(i, 0);
    }
    bucketsort(x, index, label_max);
    // add kernel values
    vector<int> count(n);
    set<int> count_index;
    int k_value;
    for (int i = 0; i < v_all; i++) {
        count_index.insert(graph_index[index_org[index[i]]]);
        count[graph_index[index_org[index[i]]]]++;
        if (i == v_all - 1 ||
                label_list(index[i], 0) != label_list(index[i + 1], 0)) {
            for (auto itr = count_index.begin(), end = count_index.end();
                    itr != end; ++itr) {
                for (auto itr2 = itr, end2 = count_index.end(); itr2 != end2;
                        ++itr2) {
                    k_value = count[*itr] * count[*itr2];
                    K_mat(*itr, *itr2) += k_value;
                    K_mat(*itr2, *itr) += k_value;
                }
                count[*itr] = 0;
            }
            count_index.clear();
        }
    }

    int v_raised_1;
    int v_raised_2;
    for (int h = 0; h < h_max; h++) {
        nei_list.setZero();

        // first put vertex label
        nei_list.col(0) = label_list.col(h);
        // second put neighbor labels
        raise = 0;
        for (int i = 0; i < n; i++) {
            fill(counter.begin(), counter.end(), 1);
            for (int j = 0; j < num_e[i]; j++) {
                v_raised_1 = E[i](j, 0) + raise;
                v_raised_2 = E[i](j, 1) + raise;
                nei_list(v_raised_1, counter[E[i](j, 0)]) = label_list(v_raised_2, h);
                nei_list(v_raised_2, counter[E[i](j, 1)]) = label_list(v_raised_1, h);
                counter[E[i](j, 0)]++;
                counter[E[i](j, 1)]++;
            }
            raise += num_v[i];
        }

        // sort each row w.r.t. neighbors
        vector<int> y(nei_list.cols() - 1);
        for (auto i = 0; i < v_all; i++) {
            for (auto j = 1L; j < nei_list.cols(); ++j) {
                y[j - 1] = nei_list(i, j);
            }
            sort(y.begin(), y.end(), greater<>());
            for (auto j = 1L; j < nei_list.cols(); ++j) {
                nei_list(i, j) = y[j - 1];
            }
        }

        // radix sort
        for (int i = 0; i < v_all; i++) {
            index[i] = i;
            index_org[i] = i;
        }
        for (auto k = nei_list.cols() - 1; k >= 0; k--) {
            for (int i = 0; i < v_all; i++) {
                x[i] = nei_list(i, k);
            }
            bucketsort(x, index, label_max);
        }

        // re-labeling and increment kernel values
        label_max++;
        for (int i = 0; i < v_all; i++) {
            label_list(index_org[index[i]], h + 1) = label_max;
            count_index.insert(graph_index[index_org[index[i]]]);
            count[graph_index[index_org[index[i]]]]++;
            if (i == v_all - 1 ||
                    (nei_list.row(index[i]) - nei_list.row(index[i + 1]))
                    .array()
                    .abs()
                    .sum() != 0) {
                for (auto itr = count_index.begin(), end = count_index.end();
                        itr != end; ++itr) {
                    for (auto itr2 = itr, end2 = count_index.end(); itr2 != end2;
                            ++itr2) {
                        k_value = count[*itr] * count[*itr2];
                        K_mat(*itr, *itr2) += k_value;
                        K_mat(*itr2, *itr) += k_value;
                    }
                    count[*itr] = 0;
                }
                count_index.clear();
                label_max++;
            }
        }
    }
    for (int i = 0; i < n; ++i) {
        K_mat(i, i) /= 2;
    }

    return K_mat;
}
