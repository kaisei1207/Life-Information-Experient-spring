#include <iostream>   
#include <fstream>    
#include <vector>     
#include <string>     
#include <cmath>      
#include <random>     
#include <algorithm>  
using namespace std;  

// ===== 塩基 → インデックス変換 =====
// A,C,G,T をそれぞれ 0,1,2,3 に変換
int base_index(char c) {
    if (c=='A') return 0; 
    if (c=='C') return 1;  
    if (c=='G') return 2;  
    if (c=='T') return 3;  
    return -1; // 不正文字
}

// ===== ランダム塩基生成 =====
// bg: バックグラウンド確率（A,C,G,Tの順）
// gen: 乱数生成器
char random_base(const vector<double>& bg, mt19937& gen) {
    uniform_real_distribution<> dis(0.0, 1.0); // 一様乱数 [0,1)
    double r = dis(gen); // 生成された乱数

    // 累積確率で塩基を決定
    if (r < bg[0]) return 'A';
    else if (r < bg[0] + bg[1]) return 'C';
    else if (r < bg[0] + bg[1] + bg[2]) return 'G';
    else return 'T';
}

// ===== ランダムDNA配列生成 =====
// length: 配列長
// bg: バックグラウンドの出現確率
string generate_random_sequence(int length, const vector<double>& bg, mt19937& gen) {
    string seq = ""; // 生成されるランダム配列
    for (int i = 0; i < length; i++) {
        seq += random_base(bg, gen);
    }
    return seq;
}

int main() {

    // ===== p値入力 =====
    double p_value; 
    cout << "p値を入力してください: ";
    cin >> p_value;

    if (p_value <= 0.0 || p_value >= 1.0) {
        cerr << "p値は 0 < p < 1 の範囲で入力してください\n";
        return 1;
    }

    // ===== 転写因子ファイル名 =====
    vector<string> filenames; // 複数の転写因子ファイル名
    string fname;             // 入力用一時変数

    cout << "転写因子ファイル名（endで終了）:" << endl;
    while (cin >> fname && fname != "end") {
        filenames.push_back(fname);
    }

    // ===== バックグラウンド出現確率 =====
    double A = 7519429; // Aの総出現回数
    double T = 7519429; // Tの総出現回数
    double C = 4637676; // Cの総出現回数
    double G = 4637676; // Gの総出現回数

    double total_bg = A + T + C + G; // 全塩基数

    // bg[b]: 塩基bのバックグラウンド出現確率
    vector<double> bg = {
        A / total_bg,
        C / total_bg,
        G / total_bg,
        T / total_bg
    };

    // ===== プロモーター配列読み込み =====
    ifstream pist("promoters"); // 入力ファイル
    if (!pist) {
        cerr << "promoters open error\n";
        return 1;
    }

    vector<string> names;      // 配列名
    vector<string> promoters;  // DNA配列本体

    string line;  // 1行読み込み
    string seq = "";   // 現在の配列
    string name = "";  // 現在の名前

    while (getline(pist, line)) {
        if (line.empty()) continue;

        if (line[0] == '>') {
            if (!seq.empty()) {
                promoters.push_back(seq); // 配列保存
                names.push_back(name);    // 名前保存
                seq = "";
            }
            name = line;
        } else {
            seq += line; // 配列連結
        }
    }

    if (!seq.empty()) {
        promoters.push_back(seq);
        names.push_back(name);
    }

    // ===== 乱数生成器 =====
    mt19937 gen(1207); // シード固定（再現性あり）

    // ===== 転写因子ごとの処理 =====
    for (string &file : filenames) {

        ifstream ist(file);
        if (!ist) {
            cerr << file << " open error\n";
            continue;
        }

        // ===== 転写因子配列 =====
        vector<string> seqs; // 結合配列の集合
        string s;

        while (getline(ist, s)) {
            if (!s.empty() && s[0] != '>') {
                seqs.push_back(s);
            }
        }

        if (seqs.empty()) continue;

        int L = seqs[0].size(); // 転写因子配列の長さ

        // ===== カウント行列 =====
        // count[b][i]: 塩基bが位置iに出た回数
        vector<vector<int>> count(4, vector<int>(L, 1)); // 疑似頻度

        for (string &seq_tf : seqs) { // 各転写因子配列
            for (int i = 0; i < L; i++) {
                switch (seq_tf[i]) {
                    case 'A': count[0][i]++; break;
                    case 'C': count[1][i]++; break;
                    case 'G': count[2][i]++; break;
                    case 'T': count[3][i]++; break;
                }
            }
        }

        // ===== 対数オッズスコア count[b][i] =====
        vector<vector<double>> pwm(4, vector<double>(L));

        for (int b = 0; b < 4; b++) {
            for (int i = 0; i < L; i++) {

                double total = count[0][i] + count[1][i] 
                             + count[2][i] + count[3][i];

                double p = (double)count[b][i] / total;

                pwm[b][i] = log(p / bg[b]); // logオッズ
            }
        }

        // ===== ランダムスコア分布 =====
        vector<double> random_scores; // スコア分布
        int num_random = 10000; // ランダム配列数

        for (int n = 0; n < num_random; n++) {

            string rand_seq = generate_random_sequence(promoters[0].size(), bg, gen);
            // ランダムDNA配列

            for (int i = 0; i <= (int)rand_seq.size() - L; i++) {

                double score = 0.0; // 部分配列のスコア

                for (int j = 0; j < L; j++) {
                    int idx = base_index(rand_seq[i + j]);
                    if (idx != -1) {
                        score += pwm[idx][j];
                    }
                }

                random_scores.push_back(score);
            }
        }

        // ===== 閾値決定 =====
        sort(random_scores.begin(), random_scores.end());

        int index = (int)((1.0 - p_value) * random_scores.size());
        // 上位(1-p)%の位置

        if (index >= random_scores.size()) index = random_scores.size() - 1;
        if (index < 0) index = 0;

        double threshold = random_scores[index]; // 閾値

        cout << "\n=== 転写因子: " << file << " ===\n";
        cout << "閾値 = " << threshold << endl;

        // ===== 実データ探索 =====
        for (int p_idx = 0; p_idx < promoters.size(); p_idx++) {

            string &sequence = promoters[p_idx]; // 対象配列

            cout << "\n" << names[p_idx] << endl;

            if (sequence.size() < (size_t)L) {
                cout << "sequence too short\n";
                continue;
            }

            // ===== スライディングウィンドウ =====
            for (int i = 0; i <= (int)sequence.size() - L; i++) {

                double score = 0.0;

                for (int j = 0; j < L; j++) {
                    int idx = base_index(sequence[i + j]);
                    if (idx != -1) {
                        score += pwm[idx][j];
                    }
                }

                // 閾値以上 → ヒット
                if (score >= threshold) {

                    string hit_seq = sequence.substr(i, L);
                    // ヒットした部分配列

                    cout << "[HIT] pos=" << i 
                         << " score=" << score 
                         << " seq=" << hit_seq << endl;
                }
            }
        }
    }
    return 0;
}