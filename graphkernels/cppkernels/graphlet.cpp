/* Copyright (c) 2018-2019 the `graphkernels` developers
 * All rights reserved.
 */

#include "graphlet.h"

#include <algorithm>
#include <numeric>

using std::accumulate;
using std::iota;
using std::max;
using std::min;
using std::vector;

using Eigen::MatrixXd;
using Eigen::VectorXd;

// ===== graphlet kernel for k = 4 ===== //
int find_min(int a, int b, int c) {
    const auto minabc = min(a, min(b, c));
    if (minabc == a) {
        if (minabc == b) {
            return (minabc == c) ? 7 : 4;
        }
        return (minabc == c) ? 5 : 1;
    }

    if (minabc == b) {
        return (minabc == c) ? 6 : 2;
    }

    return 3;
}

void card_ThreeInter(
        const vector<int>& L1,
        const vector<int>& L2,
        const vector<int>& L3,
        vector<int>& card) {
    card.resize(7);
    fill(card.begin(), card.end(), 0);
    auto i = 0;
    auto j = 0;
    auto k = 0;

    while (i < L1.size() && j < L2.size() && k < L3.size()) {
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
            default:
                {}  // FIXME: THIS SHOULD NEVER HAPPEN.
        }
    }

    if (i < L1.size() || j < L2.size() || k < L3.size()) {
        if (i >= L1.size() && j >= L2.size()) {
            card[2] += L3.size();
            card[2] -= k;
            k = L3.size();
        } else {
            if (i >= L1.size() && k >= L3.size()) {
                card[1] += L2.size();
                card[1] -= j;
                j = L2.size();
            } else {
                if (j >= L2.size() && k >= L3.size()) {
                    card[0] += L1.size();
                    card[0] -= i;
                    i = L1.size();
                } else {
                    if (i >= L1.size()) {
                        while (j < L2.size() && k < L3.size()) {
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
                        if (j >= L2.size()) {
                            while (i < L1.size() && k < L3.size()) {
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
                            if (k >= L3.size()) {
                                while (i < L1.size() && j < L2.size()) {
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
    if (i < L1.size() || j < L2.size() || k < L3.size()) {
        if (i >= L1.size() && j >= L2.size()) {
            card[2] += L3.size();
            card[2] -= k;
        } else if (i >= L1.size() && k >= L3.size()) {
            card[1] += L2.size();
            card[1] -= j;
        } else if (j >= L2.size() && k >= L3.size()) {
            card[0] += L1.size();
            card[0] -= i;
        }
    }
}

void getIndices(
        const vector<int>& o_set1,
        const vector<int>& o_set2,
        vector<int>& inter,
        vector<int>& diff1,
        vector<int>& diff2) {
    vector<int> inter_(min(o_set1.size(), o_set2.size()), -1);
    vector<int> diff1_(max(o_set1.size(), o_set2.size()), -1);
    vector<int> diff2_(max(o_set1.size(), o_set2.size()), -1);

    auto i = 0;
    auto j = 0;
    while (i < o_set1.size() && j < o_set2.size()) {
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

    if (i < o_set1.size()) {
        for (auto k = i; k < o_set1.size(); ++k) {
            diff1_[k] = o_set1[k];
        }
    } else if (j < o_set2.size()) {
        for (auto k = j; k < o_set2.size(); ++k) {
            diff2_[k] = o_set2[k];
        }
    }

    inter.clear();
    for (auto&& x : inter_) {
        if (x >= 0) {
            inter.push_back(x);
        }
    }
    diff1.clear();
    for (auto&& x : diff1_) {
        if (x >= 0) {
            diff1.push_back(x);
        }
    }
    diff2.clear();
    for (auto&& x : diff2_) {
        if (x >= 0) {
            diff2.push_back(x);
        }
    }
}

VectorXd countGraphletsFour(vector<vector<int>>& al, int freq_size) {
    VectorXd count_gr = VectorXd::Zero(freq_size);
    const auto n = al.size();
    vector<double> w = {
        1.0 / 12.0, 1.0 / 10.0, 1.0 / 8.0, 1.0 / 6.0,
        1.0 / 8.0,  1.0 / 6.0,  1.0 / 6.0, 1.0 / 4.0,
        1.0 / 4.0,  1.0 / 2.0,  0};
    vector<int> inter;
    vector<int> diff1;
    vector<int> diff2;
    vector<int> card;
    vector<double> inter_count(11);
    vector<int> v;

    double m = 0.0;
    for (auto&& vec : al) {
        m += vec.size();
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
                inter_count[0] += 0.5 * card[6];
                inter_count[1] += 0.5 * (card[3] - 1.0);
                inter_count[1] += 0.5 * (card[4] - 1.0);
                inter_count[1] += 0.5 * (card[5] - 1.0);
                inter_count[2] += 0.5 * card[0];
                inter_count[2] += 0.5 * card[1];
                inter_count[2] += card[2];
                inter_count[6] += n - accumulate(
                        card.cbegin(), card.cend(), 0.0);
                K += 0.5 * card[6] + 0.5 * (card[4] - 1.0) +
                    0.5 * (card[5] - 1.0) + card[2];
            }
            v.clear();
            v.resize(diff1.size());
            sort(diff1.begin(), diff1.end());
            sort(al[i].begin(), al[i].end());
            auto it = set_difference(diff1.cbegin(), diff1.cend(),
                    al[i].cbegin(), al[i].cend(), v.begin());
            v.resize(it - v.cbegin());
            for (auto&& k : v) {
                card_ThreeInter(al[i], al[j], al[k], card);
                inter_count[1] += 0.5 * card[6];
                inter_count[2] += 0.5 * card[3];
                inter_count[2] += 0.5 * card[4];
                inter_count[4] += 0.5 * (card[5] - 1.0);
                inter_count[3] += 0.5 * (card[0] - 2.0);
                inter_count[5] += 0.5 * card[1];
                inter_count[5] += card[2];
                inter_count[7] += n - accumulate(
                        card.cbegin(), card.cend(), 0.0);
                K += 0.5 * card[6] + 0.5 * card[4] +
                    0.5 * (card[5] - 1.0) + card[2];
            }
            v.clear();
            v.resize(diff2.size());
            sort(diff2.begin(), diff2.end());
            it = set_difference(diff2.cbegin(), diff2.cend(), v1.cbegin(),
                    v1.cend(), v.begin());
            v.resize(it - v.cbegin());
            for (auto&& k : v) {
                card_ThreeInter(al[i], al[j], al[k], card);
                inter_count[1] += 0.5 * card[6];
                inter_count[2] += 0.5 * card[3];
                inter_count[4] += 0.5 * (card[4] - 1.0);
                inter_count[2] += 0.5 * card[5];
                inter_count[5] += 0.5 * card[0];
                inter_count[3] += 0.5 * (card[1] - 2.0);
                inter_count[5] += card[2];
                inter_count[7] += n - accumulate(
                        card.cbegin(), card.cend(), 0.0);
                K += 0.5 * card[6] + 0.5 * (card[4] - 1.0) +
                    0.5 * card[5] + card[2];
            }
            inter_count[8] += m + 1.0 - v1.size() - al[i].size() - K;
            inter_count[9] +=
                (n - inter.size() - diff1.size() - diff2.size()) *
                (n - 1.0 - inter.size() - diff1.size() - diff2.size()) / 2 -
                (m + 1.0 - v1.size() - al[i].size() - K);

            for (auto k = 0L; k < count_gr.size(); ++k) {
                count_gr(k) += inter_count[k] * w[k];
            }
        }
    }

    count_gr(10) = n * (n - 1.0) * (n - 2.0) * (n - 3.0) / (4.0 * 3.0 * 2.0) -
        count_gr.head(10).sum();

    const auto csum = count_gr.sum();
    return (csum == 0.0) ? count_gr : count_gr / csum;
}

// ===== graphlet kernel for k = 3 ===== //
void getCardinality(
        const vector<int>& o_set1,
        const vector<int>& o_set2,
        vector<double>& card) {
    card.resize(3);
    fill(card.begin(), card.end(), 0.0);
    auto i = 0;
    auto j = 0;
    while (i < o_set1.size() && j < o_set2.size()) {
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
    card[0] += o_set1.size();
    card[0] -= i;

    card[1] += o_set2.size();
    card[1] -= j;
}

VectorXd countGraphletsThree(const vector<vector<int>>& al, int freq_size) {
    VectorXd count_gr = VectorXd::Zero(freq_size);
    const auto n = al.size();
    vector<double> w = {1.0 / 6.0, 1.0 / 4.0, 1.0 / 2.0};
    vector<double> card(3);

    for (auto i = 0; i < n; ++i) {
        for (auto&& j : al[i]) {
            getCardinality(al[i], al[j], card);
            count_gr(0) += w[0] * card[2];
            count_gr(1) += w[1] * (card[0] + card[1] - 2.0);
            count_gr(2) += w[2] * (
                    n - accumulate(card.cbegin(), card.cend(), 0.0));
        }
    }
    count_gr(3) = n * (n - 1.0) * (n - 2.0) / 6.0 -
        (count_gr(0) + count_gr(1) + count_gr(2));

    const auto csum = count_gr.sum();
    return (csum == 0.0) ? count_gr : count_gr / csum;
}

MatrixXd CalculateGraphletKernelThreePy(
        const vector<vector<vector<int>>>& graph_adjlist_all) {
    constexpr auto freq_size = 4;
    MatrixXd freq(graph_adjlist_all.size(), freq_size);

    for (auto i = 0; i < graph_adjlist_all.size(); ++i) {
        freq.row(i) = countGraphletsThree(graph_adjlist_all[i], freq_size);
    }
    return freq * freq.transpose();
}

MatrixXd CalculateGraphletKernelFourPy(
        vector<vector<vector<int>>>& graph_adjlist_all) {
    constexpr auto freq_size = 11;
    MatrixXd freq(graph_adjlist_all.size(), freq_size);

    for (auto i = 0; i < graph_adjlist_all.size(); ++i) {
        freq.row(i) = countGraphletsFour(graph_adjlist_all[i], freq_size);
    }
    return freq * freq.transpose();
}

MatrixXd CalculateGraphletKernelPy(
        vector<vector<vector<int>>>& graph_adjlist_all,
        int k) {
    if (k == 3) {
        return CalculateGraphletKernelThreePy(graph_adjlist_all);
    }
    return CalculateGraphletKernelFourPy(graph_adjlist_all);
}
