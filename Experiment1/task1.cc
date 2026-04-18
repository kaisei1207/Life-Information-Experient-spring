#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
using namespace std;

int main() {
    // 複数ファイル名を格納するベクタ
    vector<string> filenames;
    string fname;

    // ユーザーからファイル名を入力（"end"で終了）
    cout << "ファイル名（endで終了）:" << endl;
    while (cin >> fname && fname != "end") {
        filenames.push_back(fname);
    }

    // ===== バックグラウンド確率 q_i(x) の設定 =====
    // 各塩基の出現回数（事前に与えられた値）
    double A = 7519429, T = 7519429, C = 4637676, G = 4637676;
    double total_bg = A + T + C + G; // 全体の総数

    // 塩基ごとの確率 q（A, C, G, T の順）
    vector<double> bg(4);
    bg[0] = A / total_bg; // A の確率
    bg[1] = C / total_bg; // C の確率
    bg[2] = G / total_bg; // G の確率
    bg[3] = T / total_bg; // T の確率

    // 塩基を表示するための配列
    char base[4] = {'A','C','G','T'};

    // ===== 各ファイルについて処理 =====
    for (string &file : filenames) {

        // ファイルを開く
        ifstream ist(file);
        if (!ist) {
            cerr << file << " open error\n"; // 開けなければエラー表示
            continue;
        }

        // 配列に配列（配列の集合）として配列データを保存
        vector<string> seqs;
        string s;

        // FASTA形式を想定（">"で始まる行はヘッダなので除外）
        while (getline(ist, s)) {
            if (!s.empty() && s[0] != '>') {
                seqs.push_back(s);
            }
        }

        // 配列が空ならスキップ
        if (seqs.empty()) continue;

        // 配列の長さ（全配列同じ長さと仮定）
        int L = seqs[0].size();

        // カウント行列（4×L）を作成
        // 初期値を1にしているのは疑似頻度
        vector<vector<int>> count(4, vector<int>(L, 1));

        // ===== 塩基のカウント =====
        for (int n = 0; n < seqs.size(); n++) {
            string &seq = seqs[n];
            for (int i = 0; i < L; i++) {
                // 各位置ごとに塩基をカウント
                if (seq[i]=='A') count[0][i]++;
                if (seq[i]=='C') count[1][i]++;
                if (seq[i]=='G') count[2][i]++;
                if (seq[i]=='T') count[3][i]++;
            }
        }

        // ファイルごとの結果表示
        cout << "=== " << file << " オッズスコア log(p/q) ===\n";

        // ===== PWM（Position Weight Matrix）の計算 =====
        for (int b = 0; b < 4; b++) {
            cout << base[b] << " "; // 塩基名を表示

            for (int i = 0; i < L; i++) {
                // その位置での総カウント
                double total = count[0][i] + count[1][i]
                             + count[2][i] + count[3][i];

                // 観測確率 p（その位置での塩基頻度）
                double p = (double)count[b][i] / total;

                // バックグラウンド確率 q
                double q = bg[b];

                // オッズスコア log(p/q)
                double score = log(p / q);

                // スコアを出力
                cout << score << " ";
            }
            cout << endl;
        }

        cout << endl;
    }

    return 0;
}