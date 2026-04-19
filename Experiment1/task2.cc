#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
using namespace std;

// ===== 塩基 → インデックス変換 =====
// 文字で表される塩基を数値（配列インデックス）に変換する関数
// A=0, C=1, G=2, T=3
int base_index(char c) {
    if (c=='A') return 0;
    if (c=='C') return 1;
    if (c=='G') return 2;
    if (c=='T') return 3;
    return -1; // A,C,G,T以外（不正文字）は無視
}

int main() {

    // ===== 転写因子ファイル名の入力 =====
    vector<string> filenames; // 入力された複数の転写因子ファイル名を格納
    string fname;             // 1つのファイル名を一時的に受け取る変数

    cout << "転写因子ファイル名（endで終了）:" << endl;
    while (cin >> fname && fname != "end") {
        filenames.push_back(fname); // ファイル名を配列に追加
    }

    // ===== バックグラウンド塩基頻度 =====
    double A = 7519429; // ゲノム全体でのAの出現回数
    double T = 7519429; // Tの出現回数
    double C = 4637676; // Cの出現回数
    double G = 4637676; // Gの出現回数

    double total_bg = A + T + C + G; // 全塩基の総数

    // bg[b]: 塩基bのバックグラウンド確率
    // b=0:A, 1:C, 2:G, 3:T
    vector<double> bg = {
        A / total_bg,
        C / total_bg,
        G / total_bg,
        T / total_bg
    };

    // ===== プロモーター配列の読み込み =====
    ifstream pist("promoters"); // プロモーターファイル入力ストリーム
    if (!pist) {
        cerr << "promoters open error\n";
        return 1;
    }

    vector<string> names;      // 各配列の名前（FASTAの>行）
    vector<string> promoters;  // 各プロモーターDNA配列

    string line;  // ファイルから1行ずつ読み込むための変数
    string seq = "";   // 現在読み込み中のDNA配列
    string name = "";  // 現在の配列名

    while (getline(pist, line)) {
        if (line.empty()) continue;

        if (line[0] == '>') {
            // 新しい配列の開始
            if (!seq.empty()) {
                promoters.push_back(seq); // 1つ前の配列を保存
                names.push_back(name);    // 対応する名前を保存
                seq = ""; // 初期化
            }
            name = line; // 新しい配列名
        } else {
            seq += line; // 配列を連結（複数行対応）
        }
    }

    // 最後の配列を追加
    if (!seq.empty()) {
        promoters.push_back(seq);
        names.push_back(name);
    }

    // ===== 各転写因子ごとの処理 =====
    for (string &file : filenames) {

        ifstream ist(file); // 転写因子ファイルの入力ストリーム
        if (!ist) {
            cerr << file << " open error\n";
            continue;
        }

        // ===== 転写因子配列の読み込み =====
        vector<string> seqs; // 転写因子の結合部位配列の集合
        string s;            // 1行分の配列

        while (getline(ist, s)) {
            if (!s.empty() && s[0] != '>') {
                seqs.push_back(s); // 配列を保存
            }
        }

        if (seqs.empty()) continue;

        int L = seqs[0].size(); // モチーフの長さ（全配列同一と仮定）

        // ===== カウント行列 =====
        // count[b][i]:
        // 塩基bが位置iに何回出現したか
        vector<vector<int>> count(4, vector<int>(L, 1)); // 疑似カウント1

        // ===== 塩基の出現回数カウント =====
        for (string &seq_tf : seqs) { // seq_tf: 1つのモチーフ配列
            for (int i = 0; i < L; i++) {
                switch (seq_tf[i]) {
                    case 'A': count[0][i]++; break;
                    case 'C': count[1][i]++; break;
                    case 'G': count[2][i]++; break;
                    case 'T': count[3][i]++; break;
                }
            }
        }

        // ===== 対数オッズ行列 =====
        vector<vector<double>> pwm(4, vector<double>(L));
        // pwm[b][i]:
        // 位置iにおける塩基bの対数オッズスコア

        for (int b = 0; b < 4; b++) {
            for (int i = 0; i < L; i++) {

                double total = count[0][i] + count[1][i]
                             + count[2][i] + count[3][i];
                // 位置iにおける総カウント

                double p = (double)count[b][i] / total;
                // 観測確率（その位置での塩基頻度）

                pwm[b][i] = log(p / bg[b]);
                // 対数オッズスコア
            }
        }

        // ===== プロモーター配列の一致度 =====
        cout << "=== 転写因子: " << file << " ===\n";

        for (int p_idx = 0; p_idx < promoters.size(); p_idx++) {

            string &sequence = promoters[p_idx];
            // 現在評価しているプロモーター配列

            cout << names[p_idx] << endl;
            cout << "length = " << sequence.size() << endl;

            if (sequence.size() < (size_t)L) {
                cout << "sequence too short\n\n";
                continue;
            }

            // ===== スライディングウィンドウ =====
            for (int i = 0; i <= (int)sequence.size() - L; i++) {

                double score = 0.0;
                // 現在位置iでのマッチスコア

                for (int j = 0; j < L; j++) {
                    int idx = base_index(sequence[i + j]);
                    // 配列中の塩基をインデックスに変換

                    if (idx != -1) {
                        score += pwm[idx][j];
                        // 各位置のスコアを加算
                    }
                }

                cout << i << " " << score << ",";
                // i: 開始位置, score: 対数オッズスコアの和
            }

            cout << endl;
        }
    }

    return 0;
}