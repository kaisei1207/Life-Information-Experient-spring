#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <random>
#include <algorithm>
using namespace std;

// ===== 塩基 → インデックス =====
int base_index(char c) {
    if (c=='A') return 0;
    if (c=='C') return 1;
    if (c=='G') return 2;
    if (c=='T') return 3;
    return -1;
}

// ===== ランダム塩基 =====
char random_base(const vector<double>& bg, mt19937& gen) {
    uniform_real_distribution<> dis(0.0, 1.0);
    double r = dis(gen);

    if (r < bg[0]) return 'A';
    else if (r < bg[0] + bg[1]) return 'C';
    else if (r < bg[0] + bg[1] + bg[2]) return 'G';
    else return 'T';
}

// ===== ランダムDNA =====
string generate_random_sequence(int length, const vector<double>& bg, mt19937& gen) {
    string seq = "";
    for (int i = 0; i < length; i++) {
        seq += random_base(bg, gen);
    }
    return seq;
}

// ===== ヒット構造体 =====
struct Hit {
    int pos;
    double score;
    double p;
    string seq;
};

int main() {

    double q = 0.05;

    vector<string> filenames;
    string fname;
    cout << "転写因子ファイル名（endで終了）:" << endl;
    while (cin >> fname && fname != "end") {
        filenames.push_back(fname);
    }

    // ===== 背景確率 =====
    double A = 7519429, T = 7519429, C = 4637676, G = 4637676;
    double total_bg = A + T + C + G;

    vector<double> bg = {
        A / total_bg,
        C / total_bg,
        G / total_bg,
        T / total_bg
    };

    // ===== プロモータ読み込み =====
    ifstream pist("promoters");
    vector<string> names, promoters;

    string line, seq = "", name = "";

    while (getline(pist, line)) {
        if (line.empty()) continue;

        if (line[0] == '>') {
            if (!seq.empty()) {
                promoters.push_back(seq);
                names.push_back(name);
                seq = "";
            }
            name = line;
        } else {
            seq += line;
        }
    }

    if (!seq.empty()) {
        promoters.push_back(seq);
        names.push_back(name);
    }

    mt19937 gen(1207);

    for (string &file : filenames) {

        ifstream ist(file);
        if (!ist) continue;

        vector<string> seqs;
        string s;

        while (getline(ist, s)) {
            if (!s.empty() && s[0] != '>') {
                seqs.push_back(s);
            }
        }

        if (seqs.empty()) continue;

        int L = seqs[0].size();

        // ===== カウント =====
        vector<vector<int>> count(4, vector<int>(L, 1));

        for (string &seq_tf : seqs) {
            for (int i = 0; i < L; i++) {
                switch (seq_tf[i]) {
                    case 'A': count[0][i]++; break;
                    case 'C': count[1][i]++; break;
                    case 'G': count[2][i]++; break;
                    case 'T': count[3][i]++; break;
                }
            }
        }

        // ===== PWM =====
        vector<vector<double>> pwm(4, vector<double>(L));

        for (int b = 0; b < 4; b++) {
            for (int i = 0; i < L; i++) {
                double total = count[0][i] + count[1][i]
                             + count[2][i] + count[3][i];

                double p = (double)count[b][i] / total;
                pwm[b][i] = log(p / bg[b]);
            }
        }

        // ===== ランダム分布 =====
        vector<double> random_scores;
        int num_random = 5000;

        for (int n = 0; n < num_random; n++) {
            string rand_seq = generate_random_sequence(promoters[0].size(), bg, gen);

            for (int i = 0; i <= (int)rand_seq.size() - L; i++) {
                double score = 0.0;

                for (int j = 0; j < L; j++) {
                    int idx = base_index(rand_seq[i + j]);
                    if (idx != -1) {
                        score += pwm[idx][j];
                    }
                }
                random_scores.push_back(score);
            }
        }

        sort(random_scores.begin(), random_scores.end());

        cout << "\n=== 転写因子: " << file << " ===\n";

        for (int p_idx = 0; p_idx < promoters.size(); p_idx++) {

            cout << "\n" << names[p_idx] << endl;

            vector<Hit> hits;
            string &sequence = promoters[p_idx];

            for (int i = 0; i <= (int)sequence.size() - L; i++) {

                double score = 0.0;

                for (int j = 0; j < L; j++) {
                    int idx = base_index(sequence[i + j]);
                    if (idx != -1) {
                        score += pwm[idx][j];
                    }
                }

                auto it = lower_bound(random_scores.begin(), random_scores.end(), score);
                double p = (double)(random_scores.end() - it) / random_scores.size();

                hits.push_back({i, score, p, sequence.substr(i, L)});
            }

            // ===== p値でソート =====
            sort(hits.begin(), hits.end(), [](const Hit &a, const Hit &b) {
                return a.p < b.p;
            });

            int m = hits.size();
            int threshold_index = -1;

            for (int i = 0; i < m; i++) {
                if (hits[i].p <= (double)(i + 1) / m * q) {
                    threshold_index = i;
                }
            }

            // ===== ヒットのみ表示 =====
            if (threshold_index == -1) {
                cout << "";
            } else {
                for (int i = 0; i <= threshold_index; i++) {
                    cout << "[BH HIT] pos=" << hits[i].pos
                         << " score=" << hits[i].score
                         << " p=" << hits[i].p
                         << " seq=" << hits[i].seq << endl;
                }
            }
        }
    }
}