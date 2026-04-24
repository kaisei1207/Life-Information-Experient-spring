#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>   // iota（連番生成）
#include <random>    // shuffle用
using namespace std;

// ===== 定数 =====
// 特徴量の数（1つのデータに含まれる数値の個数）
#define NUM_FEATURES 53

// データ数（サンプル数）
#define NUM_SEQS 10000

// ============================
// データ読み込み関数
// ============================
// ファイルから
// ・特徴量の名前
// ・特徴量データ（dataset）
// ・ラベル（正解）
// を読み込む
void LoadSolubilityFile(
    const string& filename,                 // ファイル名
    vector<string>& feature_name,           // 特徴量の名前
    vector<vector<double>>& dataset,        // 特徴量データ
    vector<int>& labels)                    // 正解ラベル（0 or 1）
{
    ifstream ifs(filename); // ファイルを開く

    // ファイルが開けなかった場合はエラー
    if (!ifs) {
        cerr << "File open error\n";
        exit(1);
    }

    string tmp;

    // ===== ヘッダ読み込み =====
    // 1列目（IDなど）は使わないので読み飛ばす
    ifs >> tmp;

    // 特徴量の名前を読み込む（53個）
    for (int i = 0; i < NUM_FEATURES; i++) {
        ifs >> feature_name[i];
    }

    // 最後の列（labelなどの名前）を読み飛ばす
    ifs >> tmp;

    // ===== データ読み込み =====
    int row = 0; // 現在の行番号

    while (true) {

        string id;

        // IDが読めなければEOF（ファイル終端）なので終了
        if (!(ifs >> id)) break;

        // 各特徴量を読み込む
        for (int j = 0; j < NUM_FEATURES; j++) {
            ifs >> dataset[row][j];
        }

        // ラベル（0 or 1）を読み込む
        ifs >> labels[row];

        row++; // 次の行へ
    }

    // 実際に読み込んだデータ数を表示
    cout << "Loaded: " << row << endl;
}

// ============================
// データ分割関数
// ============================
// datasetを
// ・訓練データ
// ・テストデータ
// にランダムに分ける
void DivideDataset(
    const vector<vector<double>>& dataset,  // 全データ
    const vector<int>& labels,              // 全ラベル
    vector<vector<double>>& training_dataset, // 訓練データ
    vector<int>& training_labels,
    vector<vector<double>>& test_dataset,     // テストデータ
    vector<int>& test_labels,
    double test_ratio)                      // テストデータの割合
{
    int N = dataset.size(); // データ数

    // インデックス（0〜N-1）を入れる配列
    vector<int> indices(N);

    // 0,1,2,...,N-1 を作る
    // → データの「番号」を管理するため
    iota(indices.begin(), indices.end(), 0);

    // ランダムに並び替え
    // → データをランダムに分割するため
    mt19937 gen(1207); // 乱数（固定値で再現可能）
    shuffle(indices.begin(), indices.end(), gen);

    // テストデータの個数
    int test_size = (int)(N * test_ratio);

    // ===== テストデータ =====
    for (int i = 0; i < test_size; i++) {
        int idx = indices[i]; // シャッフルされた番号
        test_dataset.push_back(dataset[idx]);
        test_labels.push_back(labels[idx]);
    }

    // ===== 訓練データ =====
    for (int i = test_size; i < N; i++) {
        int idx = indices[i];
        training_dataset.push_back(dataset[idx]);
        training_labels.push_back(labels[idx]);
    }
}

// ============================
// 評価関数
// ============================
// すべて「1」と予測する単純モデルの性能を計算
void Evaluation(
    const vector<vector<double>>& X, // 特徴量（今回は未使用）
    const vector<int>& y)            // 正解ラベル
{
    // 混同行列の要素
    int TP = 0, FP = 0, FN = 0, TN = 0;

    // 全データについて評価
    for (int i = 0; i < y.size(); i++) {

        int pred = 1; // 常に「1」と予測（ベースライン）

        // 条件ごとにカウント
        if (pred == 1 && y[i] == 1) TP++; // 正しく1
        else if (pred == 1 && y[i] == 0) FP++; // 間違って1
        else if (pred == 0 && y[i] == 1) FN++; // 間違って0
        else if (pred == 0 && y[i] == 0) TN++; // 正しく0
    }

    // ===== 評価指標 =====
    double accuracy = (double)(TP + TN) / y.size(); // 正解率

    double precision = (TP + FP == 0) ? 0 :
        (double)TP / (TP + FP); // 適合率

    double recall = (TP + FN == 0) ? 0 :
        (double)TP / (TP + FN); // 再現率

    double fscore = (precision + recall == 0) ? 0 :
        2 * precision * recall / (precision + recall); // F値

    // 結果表示
    cout << "Accuracy: " << accuracy << endl;
    cout << "Precision: " << precision << endl;
    cout << "Recall: " << recall << endl;
    cout << "F-score: " << fscore << endl;

    // 混同行列
    cout << "Confusion Matrix" << endl;
    cout << "TP: " << TP << " FP: " << FP << endl;
    cout << "FN: " << FN << " TN: " << TN << endl;
}

// ============================
// main関数
// ============================
int main(void){

    // ===== 特徴量の名前 =====
    // 53個の特徴量名を格納
    vector<string> feature_name(NUM_FEATURES, "");

    // ===== データ本体 =====
    // 10000行 × 53列の行列
    // dataset[i][j]：i番目のデータのj番目の特徴量
    vector<vector<double> > dataset(NUM_SEQS, vector<double>(NUM_FEATURES, 0.0));

    // ===== ラベル =====
    // 各データの正解（0 or 1）
    vector<int> labels (NUM_SEQS);

    // ファイルからデータ読み込み
    LoadSolubilityFile("protein_solubility_dataset.txt",
        feature_name, dataset, labels);

    // ===== 分割後データ =====
    vector<vector<double> > training_dataset; // 学習用
    vector<int> training_labels;

    vector<vector<double> > test_dataset;     // 評価用
    vector<int> test_labels;

    double test_ratio = 0.2; // 20%をテストデータに

    // データ分割
    DivideDataset(dataset, labels,
        training_dataset, training_labels,
        test_dataset, test_labels,
        test_ratio);

    // ===== ベースライン評価 =====
    // 「全部1と予測するモデル」で性能確認
    Evaluation(test_dataset, test_labels);

    return 0;
}