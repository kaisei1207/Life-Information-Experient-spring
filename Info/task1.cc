#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
using namespace std;

int main() {
    vector<string> filenames;
    string fname;

    cout << "ファイル名（endで終了）:" << endl;
    while (cin >> fname && fname != "end") {
        filenames.push_back(fname);
    }

    // バックグラウンド確率 q_i(x)
    double A = 7519429, T = 7519429, C = 4637676, G = 4637676;
    double total_bg = A + T + C + G;

    vector<double> bg(4);
    bg[0] = A / total_bg; // A
    bg[1] = C / total_bg; // C
    bg[2] = G / total_bg; // G
    bg[3] = T / total_bg; // T

    char base[4] = {'A','C','G','T'};

    for (string &file : filenames) {

        ifstream ist(file);
        if (!ist) {
            cerr << file << " open error\n";
            continue;
        }

        vector<string> seqs;
        string s;
        while (getline(ist, s)) {
            if (!s.empty() && s[0] != '>') {
                seqs.push_back(s);
            }
        }

        if (seqs.empty()) continue;

        int L = seqs[0].size();
        vector<vector<int>> count(4, vector<int>(L, 1));

        // 塩基カウント
        for (int n = 0; n < seqs.size(); n++) {
            string &seq = seqs[n];
            for (int i = 0; i < L; i++) {
                if (seq[i]=='A') count[0][i]++;
                if (seq[i]=='C') count[1][i]++;
                if (seq[i]=='G') count[2][i]++;
                if (seq[i]=='T') count[3][i]++;
            }
        }

        cout << "=== " << file << " オッズスコア log(p/q) ===\n";

        // PWM計算
        for (int b = 0; b < 4; b++) {
            cout << base[b] << " ";
            for (int i = 0; i < L; i++) {
                double total = count[0][i] + count[1][i]
                             + count[2][i] + count[3][i];

                double p = (double)count[b][i] / total;
                double q = bg[b];
                double score = log(p / q);

                cout << score << " ";
            }
            cout << endl;
        }

        cout << endl;
    }

    return 0;
}