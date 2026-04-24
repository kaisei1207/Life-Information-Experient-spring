#include <iostream>   
#include <fstream>   
#include <vector>     
#include <string>     
#include <algorithm>  
#include <numeric>    
#include <random>     
using namespace std;

// ===== 定数 =====
#define NUM_FEATURES 53   // 特徴量の数
#define NUM_SEQS 10000    // データ数

// ============================
// 決定木ノード（1回の分岐）
// ============================
struct TreeNode{
    int feature_id;        // どの特徴量で分岐するか（列番号）
    double threshold;      // 分割する閾値
    int left_class_id;     // 左側（<=）の予測ラベル
    int right_class_id;    // 右側（>）の予測ラベル
};

// ============================
// データ読み込み
// ============================
void LoadSolubilityFile(
    const string& filename,                // ファイル名
    vector<string>& feature_name,          // タンパク質名
    vector<vector<double>>& dataset,       // 特徴量データ
    vector<int>& labels)                   // 可溶性ラベル
{
    ifstream ifs(filename);               
    if (!ifs) {                            
        cerr << "File open error\n";
        exit(1);
    }

    string tmp;

    ifs >> tmp; // 先頭の不要な文字（IDなど）を読み飛ばす

    // 特徴量名を読み込む
    for (int i = 0; i < NUM_FEATURES; i++) {
        ifs >> feature_name[i];
    }

    ifs >> tmp; // 最後のラベル列名を読み飛ばす

    int row = 0; // 行番号

    while (true) {
        string id;

        if (!(ifs >> id)) break; // 読めなければ終了（EOF）

        // 特徴量を読み込む
        for (int j = 0; j < NUM_FEATURES; j++) {
            ifs >> dataset[row][j];
        }

        ifs >> labels[row]; // 可溶性ラベル（0 or 1）

        row++;
    }

    cout << "Loaded: " << row << endl;
}

// ============================
// データ分割（学習用とテスト用）
// ============================
void DivideDataset(
    const vector<vector<double>>& dataset,       // dataset : 特徴量データ（全データ）
    const vector<int>& labels,                   // labels : 正解ラベル（全データ）
    vector<vector<double>>& training_dataset,    // training_dataset : 学習用データ（出力）
    vector<int>& training_labels,                // training_labels : 学習用ラベル（出力）
    vector<vector<double>>& test_dataset,        // test_dataset : テスト用データ（出力）
    vector<int>& test_labels,                    // test_labels : テスト用ラベル（出力）
    double test_ratio)                           // test_ratio : テストデータの割合（例：0.2 → 20%）
{
    int N = dataset.size(); // データの総数（行数）

    vector<int> indices(N); // インデックスを格納する配列
                            // 各データの「番号」を持つためのもの

    iota(indices.begin(), indices.end(), 0); // indices = {0,1,2,3,...,N-1} を作る
                                             // →各データの位置を表す番号
    mt19937 gen(1207); // 乱数生成器(sheed値)

    shuffle(indices.begin(), indices.end(), gen);  // indicesの順番をランダムに並び替える
                                                   // → データをランダムに分割するため
    int test_size = (int)(N * test_ratio);  // テストデータの数
                                            // 例：N=10000, ratio=0.2 → 2000
    // ===== テストデータ作成 =====
    for (int i = 0; i < test_size; i++) {
        int idx = indices[i];               // シャッフルされたインデックスから取り出す
                                            // → ランダムに選ばれたデータ番号
        test_dataset.push_back(dataset[idx]); // 対応する特徴量をテストデータへ追加
        test_labels.push_back(labels[idx]);   // 対応するラベルもテスト用へ追加
    }

    // ===== トレーニングデータ作成 =====
    for (int i = test_size; i < N; i++) {
        int idx = indices[i]; // 残りのデータ（テストに使わなかった部分）
        training_dataset.push_back(dataset[idx]); // 学習用特徴量として追加
        training_labels.push_back(labels[idx]);   // 学習用ラベルとして追加
    }
}

// ============================
// Gini不純度
// ============================
// count1 : ラベル1の個数
// total  : データ数
// → 「どれくらいラベルが混ざっているか」を表す指標
double Gini(int count1, int total) {
    if (total == 0) return 0; // データが存在しない場合は不純度0とする（安全処理）
    double p = (double)count1 / total;  // ラベル1の割合
    return 2 * p * (1 - p);  // Gini不純度の計算
}

// ============================
// 多数決
// ============================
// count1 : ラベル1の個数
// total  : データ総数
// → 「このグループは0と1どちらが多いか」で予測ラベルを決める
int Majority(int count1, int total) {
    int count0 = total - count1; // ラベル0の数
    // ラベル1の方が多い（または同数）の場合
    if (count1 >= count0) {
        return 1; // 1を予測
    } else {
        return 0; // 0を予測
    }
}

// ============================
// 決定木学習
// ============================
// X : 特徴量データ（行＝データ、列＝特徴量）
// y : 可溶性ラベル
// node : 学習結果（どの特徴量・どの閾値かを保存）
void TrainDecisionNode(
    const vector<vector<double>>& X,
    const vector<int>& y,
    TreeNode& node)
{
    int N = X.size(); // データ数（行数）
    double best_gini = 1.0; // 今までで最も良い（小さい）Giniを記録
                            // 初期値は最大(1.0)にしておく

    // ===== 全ての特徴量を試す =====
    for (int f = 0; f < NUM_FEATURES; f++) {
        vector<double> values(N); // この特徴量fの値だけを取り出す配列

        // 各データの「f番目の特徴量」を取り出す
        for (int i = 0; i < N; i++) {
            values[i] = X[i][f];
        }
        sort(values.begin(), values.end()); // 値を小さい順に並べる
                                            // → 閾値候補を取りやすくするため
        // ===== 閾値候補を100個試す =====
        for (int p = 1; p < 100; p++) {
            int idx = (int)(N * p / 100.0); // データのp%地点のインデックス
            double threshold = values[idx]; // 閾値（この値で分割）
            int L_total = 0, R_total = 0;   // 左右のデータ数
            int L_1 = 0, R_1 = 0;           // 左右でラベル1の数
            // ===== 全データを左右に振り分け =====
            for (int i = 0; i < N; i++) {
                if (X[i][f] <= threshold) {
                    // 閾値以下 → 左
                    L_total++;              // 左の総数++
                    if (y[i] == 1) L_1++;   // ラベル1ならカウント
                }
                else {
                    // 閾値より大 → 右
                    R_total++;
                    if (y[i] == 1) R_1++;
                }
            }

            // 片側にデータがない場合はスキップ
            if (L_total == 0 || R_total == 0) continue;

            // ===== Gini不純度を計算 =====
            double GL = Gini(L_1, L_total); // 左の不純度
            double GR = Gini(R_1, R_total); // 右の不純度

            // 全体の不純度（重み付き平均）
            double G =
                (double)L_total / N * GL +
                (double)R_total / N * GR;

            // ===== 最良の分割を更新 =====
            if (G < best_gini) {
                best_gini = G;        // 最小値更新
                node.feature_id = f;  // どの特徴量か
                node.threshold = threshold; // 閾値

                // 左右それぞれ多数決でラベル決定
                node.left_class_id  = Majority(L_1, L_total);
                node.right_class_id = Majority(R_1, R_total);
            }
        }
    }
}

// ============================
// 予測
// ============================
// node : 学習済みの決定木（1回だけ分岐）
// x    : 1つのデータ（特徴量ベクトル）
int Predict(const TreeNode& node, const vector<double>& x)
{
    double value = x[node.feature_id]; // 分岐に使う特徴量の値を取り出す
    if (value <= node.threshold) {
        // 閾値以下なら左へ
        return node.left_class_id;
    } else {
        // 閾値より大なら右へ
        return node.right_class_id;
    }
}

// ============================
// 評価
// ============================
// TP: 正しく1と予測
// FP: 本当は0なのに1と予測
// FN: 本当は1なのに0と予測
// TN: 正しく0と予測
void Evaluation(
    const TreeNode& node,
    const vector<vector<double>>& X,
    const vector<int>& y)
{
    int TP = 0, FP = 0, FN = 0, TN = 0; // 混同行列の各要素を初期化

    // ===== 全データで評価 =====
    for (int i = 0; i < X.size(); i++) {

        int pred = Predict(node, X[i]); // 予測ラベル

        // ===== 結果を分類 =====
        if (pred == 1 && y[i] == 1) TP++;       // 正しく1
        else if (pred == 1 && y[i] == 0) FP++;  // 誤って1
        else if (pred == 0 && y[i] == 1) FN++;  // 誤って0
        else TN++;                              // 正しく0
    }

    // 正解率
    double acc = (double)(TP + TN) / X.size();

    // ===== Precision（適合率）=====
    double prec;
    if (TP + FP == 0) {
        prec = 0; // 分母0対策
    } else {
        prec = (double)TP / (TP + FP);
    }

    // ===== Recall（再現率）=====
    double rec;
    if (TP + FN == 0) {
        rec = 0;
    } else {
        rec = (double)TP / (TP + FN);
    }

    // ===== F1スコア =====
    double f1;
    if (prec + rec == 0) {
        f1 = 0;
    } else {
        f1 = 2 * prec * rec / (prec + rec);
    }

    // ===== 結果表示 =====
    cout << "Accuracy: " << acc << endl;
    cout << "Precision: " << prec << endl;
    cout << "Recall: " << rec << endl;
    cout << "F-score: " << f1 << endl;
    cout << "feature_id: " << node.feature_id << endl;
    cout << "threshold: " << node.threshold << endl;

    // ===== 混同行列表示 =====
    cout << endl;
    cout << "Confusion Matrix" << endl;
    cout << "TP: " << TP << " FP: " << FP << endl;
    cout << "FN: " << FN << " TN: " << TN << endl;
}

// ============================
// main関数
// ============================
int main(){

    // ===== 特徴量名を格納する配列 =====
    vector<string> feature_name(NUM_FEATURES);
    // 各特徴量の名前を保存する（例：長さ、電荷など）
    // 要素数はNUM_FEATURES（=53）

    // ===== データ本体（特徴量） =====
    vector<vector<double>> dataset(NUM_SEQS, vector<double>(NUM_FEATURES));
    // ↑ 2次元配列
    // ・行：データ数（NUM_SEQS = 10000個のタンパク質）
    // ・列：特徴量（53個）
    // dataset[i][j] → i番目データのj番目特徴量

    // ===== 正解ラベル =====
    vector<int> labels(NUM_SEQS);
    // ↑ 各データの正解（0 or 1）
    // 例：溶けない=0、溶ける=1

    // ===== ファイルからデータ読み込み =====
    LoadSolubilityFile("protein_solubility_dataset.txt",
        feature_name, dataset, labels);
    // ・feature_name（特徴量名）
    // ・dataset（特徴量データ）
    // ・labels（正解ラベル）

    // ===== 学習用・テスト用データの入れ物 =====
    vector<vector<double>> training_dataset, test_dataset; // 特徴量データ（学習用 / テスト用）
    vector<int> training_labels, test_labels; // ラベル（学習用 / テスト用）

    // ===== データ分割 =====
    DivideDataset(dataset, labels,
        training_dataset, training_labels,  // 80% → training（学習用）
        test_dataset, test_labels, 0.2);    // 20% → test（テスト用）

    // ===== 決定木（深さ１）を作る =====
    TreeNode decision_tree; // 「どの特徴量・どの閾値で分けるか」を保存する箱

    // ===== 学習 =====
    TrainDecisionNode(training_dataset, training_labels, decision_tree);
    // 学習データを使って「最も良い分割（特徴量と閾値）」を探し、decision_tree に結果を保存する

    // ===== 評価 =====
    Evaluation(decision_tree, test_dataset, test_labels);
    // テストデータを使って性能を測る
    // Accuracy（正解率）
    // Precision（適合率）
    // Recall（再現率）
    // F1スコア

    return 0;
}