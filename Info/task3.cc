#include <iostream>   
#include <fstream>    
#include <vector>     
#include <string>     
#include <cmath>      
#include <random>     
#include <algorithm>  
using namespace std;  

// 塩基（A,C,G,T）を数値（0,1,2,3）に変換する関数
int base_index(char c) {
    if (c=='A') return 0; 
    if (c=='C') return 1;  
    if (c=='G') return 2;  
    if (c=='T') return 3;  
    return -1;             
}

// 背景確率 bg に従ってランダムな塩基を1つ生成する関数
char random_base(const vector<double>& bg, mt19937& gen) {
    uniform_real_distribution<> dis(0.0, 1.0); // 0〜1の一様乱数
    double r = dis(gen); // 乱数を生成

    // 累積確率で塩基を決定
    if (r < bg[0]) return 'A';                         // A
    else if (r < bg[0] + bg[1]) return 'C';            // C
    else if (r < bg[0] + bg[1] + bg[2]) return 'G';    // G
    else return 'T';                                   // T
}

// ランダムDNA配列を生成
string generate_random_sequence(int length, const vector<double>& bg, mt19937& gen) {
    string seq = "";
    for (int i = 0; i < length; i++) {
        seq += random_base(bg, gen);
    }
    return seq;
}

int main() {

    // p値入力
    double p_value;
    cout << "p値を入力してください: ";
    cin >> p_value;

    if (p_value <= 0.0 || p_value >= 1.0) {
        cerr << "p値は 0 < p < 1 の範囲で入力してください\n";
        return 1;
    }

    // 転写因子ファイル名入力
    vector<string> filenames;
    string fname;

    cout << "転写因子ファイル名（endで終了）:" << endl;
    while (cin >> fname && fname != "end") {
        filenames.push_back(fname);
    }

    // 背景確率
    double A = 7519429, T = 7519429, C = 4637676, G = 4637676;
    double total_bg = A + T + C + G;

    vector<double> bg = {
        A / total_bg,
        C / total_bg,
        G / total_bg,
        T / total_bg
    };

    // プロモーター読み込み
    ifstream pist("promoters");
    if (!pist) {
        cerr << "promoters open error\n";
        return 1;
    }

    vector<string> names;
    vector<string> promoters;

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

    // 乱数生成器
    mt19937 gen(1207);

    // 転写因子ごとの処理
    for (string &file : filenames) {

        ifstream ist(file);
        if (!ist) {
            cerr << file << " open error\n";
            continue;
        }

        // モチーフ配列読み込み
        vector<string> seqs;
        string s;

        while (getline(ist, s)) {
            if (!s.empty() && s[0] != '>') {
                seqs.push_back(s);
            }
        }

        if (seqs.empty()) continue;

        int L = seqs[0].size();

        // カウント行列
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

        // PWM 
        vector<vector<double>> pwm(4, vector<double>(L));

        for (int b = 0; b < 4; b++) {
            for (int i = 0; i < L; i++) {
                double total = count[0][i] + count[1][i] 
                             + count[2][i] + count[3][i];

                double p = (double)count[b][i] / total;

                pwm[b][i] = log(p / bg[b]); // 対数オッズ
            }
        }

        // ランダムスコア分布
        vector<double> random_scores;
        int num_random = 5000; //ランダムDNA配列の本数

        for (int n = 0; n < num_random; n++) {

            // ランダムDNA配列１本生成
            string rand_seq = generate_random_sequence(promoters[0].size(), bg, gen);

            // 部分配列を1文字ずつずらす
            for (int i = 0; i <= (int)rand_seq.size() - L; i++) {

                double score = 0.0;

                for (int j = 0; j < L; j++) {
                    int idx = base_index(rand_seq[i + j]);
                    if (idx != -1) {
                        score += pwm[idx][j]; // PWMスコアを加算
                    }
                }

                random_scores.push_back(score); // 分布を作るための配列
            }
        }

        // 閾値決定
        sort(random_scores.begin(), random_scores.end());

        // p値に基づいて閾値の位置を決定
        int index = (int)((1.0 - p_value) * random_scores.size());

        if (index >= random_scores.size()) index = random_scores.size() - 1;
        if (index < 0) index = 0;

        double threshold = random_scores[index];

        cout << "\n=== 転写因子: " << file << " ===\n";
        cout << "閾値 = " << threshold << endl;

        // 実データ探索
        for (int p_idx = 0; p_idx < promoters.size(); p_idx++) {

            string &sequence = promoters[p_idx];

            cout << "\n" << names[p_idx] << endl;
            cout << "length = " << sequence.size() << endl;

            if (sequence.size() < (size_t)L) {
                cout << "sequence too short\n";
                continue;
            }

            // スライディングウィンドウで探索
            // 長さLの部分配列を全部調べる
            for (int i = 0; i <= (int)sequence.size() - L; i++) {

                double score = 0.0;

                for (int j = 0; j < L; j++) {
                    int idx = base_index(sequence[i + j]);
                    if (idx != -1) {
                        score += pwm[idx][j];
                    }
                }

                // 閾値以上ならヒット
                //「ランダムでは出にくいスコア」= 意味のある一致
                if (score >= threshold) {
                    // ヒットした部分配列（長さL）を取得
                    string hit_seq = sequence.substr(i, L);
                    // sequence：プロモーター配列
                    // i：開始位置
                    // L：モチーフ長
                    
                    cout << "[HIT] pos=" << i 
                    << " score=" << score 
                    << " seq=" << hit_seq << endl;
                }
            }
        }
    }
    return 0;
}