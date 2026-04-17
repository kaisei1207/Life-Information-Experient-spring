#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
using namespace std;

// ===== 塩基 → インデックス変換 =====
// PWMやカウント配列で使うため、塩基を数値に対応させる
// A=0, C=1, G=2, T=3
int base_index(char c) {
    if (c=='A') return 0;
    if (c=='C') return 1;
    if (c=='G') return 2;
    if (c=='T') return 3;
    return -1; // A,C,G,T以外（無視する）
}

int main() {

    // ===== 転写因子ファイル名の入力 =====
    // 複数の転写因子（モチーフ）を処理するため、ファイル名を配列に保存
    vector<string> filenames;
    string fname;  // ← 1回分の入力を一時的に受け取る変数

    cout << "転写因子ファイル名（endで終了）:" << endl;
    while (cin >> fname && fname != "end") {
        filenames.push_back(fname); // 入力された名前を保存
    }

    // ===== バックグラウンド塩基頻度 =====
    // ゲノム全体での塩基の出現頻度（事前に与えられている）
    double A = 7519429, T = 7519429, C = 4637676, G = 4637676;
    double total_bg = A + T + C + G;

    // 各塩基の確率 P(bg)
    vector<double> bg = {
        A / total_bg,
        C / total_bg,
        G / total_bg,
        T / total_bg
    };

    // ===== プロモーター配列の読み込み =====
    // FASTA形式：
    // >名前
    // 配列
    ifstream pist("promoters");
    if (!pist) {
        cerr << "promoters open error\n";
        return 1;
    }

    vector<string> names;      // 配列名（>xxx）
    vector<string> promoters;  // 配列本体

    string line, seq = "", name = "";

    while (getline(pist, line)) {
        if (line.empty()) continue;

        if (line[0] == '>') {
            // 新しい配列の開始
            if (!seq.empty()) {
                promoters.push_back(seq); // これまでの配列を保存
                names.push_back(name);
                seq = "";
            }
            name = line; // 配列名
        } else {
            seq += line; // 配列を連結
        }
    }

    // 最後の配列も忘れずに追加
    if (!seq.empty()) {
        promoters.push_back(seq);
        names.push_back(name);
    }

    // ===== 各転写因子ごとに処理 =====
    for (string &file : filenames) {

        ifstream ist(file);
        if (!ist) {
            cerr << file << " open error\n";
            continue;
        }

        // ===== 転写因子の結合配列（モチーフ）を読み込み =====
        vector<string> seqs;
        string s;
        while (getline(ist, s)) {
            if (!s.empty() && s[0] != '>') {
                seqs.push_back(s);
            }
        }

        if (seqs.empty()) continue;

        int L = seqs[0].size(); // モチーフの長さ

        // ===== カウント行列 =====
        // count[塩基][位置]
        // 初期値1は「疑似カウント（0割防止＆安定化）」
        vector<vector<int>> count(4, vector<int>(L, 1));

        // ===== 各位置ごとの塩基出現数をカウント =====
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

        // ===== PWM（Position Weight Matrix）作成 =====
        // 各位置で「その塩基がどれくらい出やすいか」を表す
        vector<vector<double>> pwm(4, vector<double>(L));

        for (int b = 0; b < 4; b++) {
            for (int i = 0; i < L; i++) {

                // その位置の総カウント
                double total = count[0][i] + count[1][i]
                             + count[2][i] + count[3][i];

                // 観測確率 P(x)
                double p = (double)count[b][i] / total;

                // ログオッズ：
                // 「背景よりどれだけ出やすいか」
                pwm[b][i] = log(p / bg[b]);
            }
        }

        // ===== プロモーター配列に対してスコア計算 =====
        cout << "=== 転写因子: " << file << " ===\n";

        for (int p_idx = 0; p_idx < promoters.size(); p_idx++) {
            string &sequence = promoters[p_idx];

            cout << names[p_idx] << endl;
            cout << "length = " << sequence.size() << endl;

            // モチーフより短い場合はスキップ
            if (sequence.size() < (size_t)L) {
                cout << "sequence too short\n\n";
                continue;
            }

            // ===== スライディングウィンドウ =====
            // 長さLの窓を1文字ずつずらしてスコアを計算
            for (int i = 0; i <= (int)sequence.size() - L; i++) {

                double score = 0.0;

                // 各位置のPWMを足し合わせる
                for (int j = 0; j < L; j++) {
                    int idx = base_index(sequence[i + j]);

                    if (idx != -1) {
                        score += pwm[idx][j];
                    }
                }

                // i: 開始位置, score: マッチ度
                cout << i << " " << score << ",";
            }

            cout << endl;
        }
    }

    return 0;
}