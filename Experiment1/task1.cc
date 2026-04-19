#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
using namespace std;

int main() {
    // 入力されたファイル名を格納する配列
    vector<string> filenames;
    string fname; // 1つずつ読み込むファイル名

    // ユーザーからファイル名を入力（"end"で終了）
    cout << "ファイル名（endで終了）:" << endl;
    while (cin >> fname && fname != "end") {
        filenames.push_back(fname); // 入力されたファイル名を保存
    }

    // ===== バックグラウンド確率 q_i(x) の設定 =====
    // 各塩基の出現回数（事前に与えられた全体統計）
    double A = 7519429; // 塩基Aの総出現回数
    double T = 7519429; // 塩基Tの総出現回数
    double C = 4637676; // 塩基Cの総出現回数
    double G = 4637676; // 塩基Gの総出現回数

    double total_bg = A + T + C + G; // 全塩基の総数

    // 塩基ごとのバックグラウンド確率 q
    // インデックス: 0=A, 1=C, 2=G, 3=T
    vector<double> bg(4);
    bg[0] = A / total_bg; // q(A)
    bg[1] = C / total_bg; // q(C)
    bg[2] = G / total_bg; // q(G)
    bg[3] = T / total_bg; // q(T)

    // 表示用：インデックスと塩基の対応
    char base[4] = {'A','C','G','T'};

    // ===== 各ファイルについて処理 =====
    for (string &file : filenames) {

        ifstream ist(file); // 入力ファイルストリーム
        if (!ist) {
            cerr << file << " open error\n"; // 開けない場合
            continue;
        }

        // DNA配列を格納する配列
        vector<string> seqs;
        string s; // 1行分の配列

        // 1行1配列として読み込む（空行は無視）
        while (getline(ist, s)) {
            if (!s.empty()) {
                seqs.push_back(s); // 配列を保存
            }
        }

        // データが無ければ次へ
        if (seqs.empty()) continue;

        // 配列長 L（すべて同じ長さと仮定）
        int L = seqs[0].size();

        // カウント行列（4×L）
        // count[b][i] = 塩基bが位置iに何回出現したか
        // 初期値1は疑似頻度
        vector<vector<int>> count(4, vector<int>(L, 1));

        // ===== 塩基のカウント =====
        for (int n = 0; n < seqs.size(); n++) {
            string &seq = seqs[n]; // n番目の配列

            for (int i = 0; i < L; i++) {
                // 各位置 i における塩基をカウント
                if (seq[i]=='A') count[0][i]++;
                if (seq[i]=='C') count[1][i]++;
                if (seq[i]=='G') count[2][i]++;
                if (seq[i]=='T') count[3][i]++;
            }
        }

        // 結果表示
        cout << "=== " << file << " オッズスコア log(p/q) ===\n";

        // ===== PWM（Position Weight Matrix）の計算 =====
        for (int b = 0; b < 4; b++) {
            cout << base[b] << " "; // 塩基名

            for (int i = 0; i < L; i++) {

                // 位置 i における総カウント数
                double total = count[0][i] + count[1][i]
                             + count[2][i] + count[3][i];

                // p = 観測確率（その位置での塩基bの頻度）
                double p = (double)count[b][i] / total;

                // q = バックグラウンド出現確率
                double q = bg[b];

                // logオッズスコア log(p/q)
                double score = log(p / q);

                // 出力
                cout << score << " ";
            }
            cout << endl;
        }

        cout << endl;
    }

    return 0;
}
