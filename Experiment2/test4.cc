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
    int feature_id;        // 使用する特徴量番号
    double threshold;      // 分岐の閾値

    // 左右へ進んだ場合の予測ラベル
    // 1 : 可溶性を示す
    // 0 : 可溶性を示さない
    int left_class_id;
    int right_class_id;
};


// ============================
// データ読み込み
// ============================
void LoadSolubilityFile(
    const string& filename,             // ファイル名
    vector<string>& feature_name,       // 特徴量名
    vector<vector<double>>& dataset,    // 特徴量データ
    vector<int>& labels)                // 正解ラベル
{
    ifstream ifs(filename);

    if (!ifs) {
        cerr << "File open error\n";
        exit(1);
    }

    string tmp;

    ifs >> tmp; // ID列をスキップ

    // 特徴量名読み込み
    for (int i = 0; i < NUM_FEATURES; i++) {
        ifs >> feature_name[i];
    }

    ifs >> tmp; // ラベル列名スキップ

    int row = 0;

    while (true) {

        string id;

        if (!(ifs >> id))
            break;

        // 特徴量
        for (int j = 0; j < NUM_FEATURES; j++) {
            ifs >> dataset[row][j];
        }

        // ラベル
        ifs >> labels[row];

        row++;
    }

    cout << "Loaded: " << row << endl;
}


// ============================
// データ分割（学習/テスト）
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

    // ===== 学習データ作成 =====
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
double Gini(int count1, int total)
{
    if (total == 0)
        return 0; // データが存在しない場合は不純度0とする（安全処理）

    double p = (double)count1 / total;  // ラベル1の割合

    return 2 * p * (1 - p);  // Gini不純度の計算
}


// ============================
// 多数決
// ============================
// count1 : ラベル1の個数
// total  : データ総数
// → 「このグループは0と1どちらが多いか」で予測ラベルを決める
int Majority(int count1, int total)
{
    int count0 = total - count1; // ラベル0の数（全体 - ラベル1）

    // ラベル1の方が多い（または同数）の場合
    if (count1 >= count0) {
        return 1; // 1を予測
    }
    else {
        return 0; // 0を予測
    }
}


// ============================
// ノード学習（深さ1）
// ============================
// X : 特徴量データ（行＝データ、列＝特徴量）
// y : 正解ラベル
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

                    if (y[i] == 1)
                        L_1++;   // ラベル1ならカウント
                }
                else {

                    // 閾値より大 → 右
                    R_total++;

                    if (y[i] == 1)
                        R_1++;
                }
            }

            // 片側にデータがない場合はスキップ
            if (L_total == 0 || R_total == 0)
                continue;

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
// 決定木（深さ2）学習
// ============================
// tree[0] : 親ノード
// tree[1] : 左側の子ノード
// tree[2] : 右側の子ノード
void TrainDecisionTree(
    const vector<vector<double>>& X,
    const vector<int>& y,
    vector<TreeNode>& tree)
{
    // ===== 親ノードの学習 =====
    // 全データを使って「最も良い分割」を探す
    TrainDecisionNode(X, y, tree[0]);

    // ===== データを左右に分割するための準備 =====
    vector<vector<double>> X_left, X_right; // 特徴量
    vector<int> y_left, y_right;            // ラベル

    // ===== 親ノードの条件でデータを分割 =====
    for (int i = 0; i < X.size(); i++) {

        // 親ノードの条件（feature_id + threshold）で分岐
        if (X[i][tree[0].feature_id] <= tree[0].threshold) {

            // 左側に入るデータ
            X_left.push_back(X[i]);
            y_left.push_back(y[i]);

        }
        else {

            // 右側に入るデータ
            X_right.push_back(X[i]);
            y_right.push_back(y[i]);
        }
    }

    // ===== 左側の子ノード学習 =====
    if (!X_left.empty())
        TrainDecisionNode(X_left, y_left, tree[1]);

    // ===== 右側の子ノード学習 =====
    if (!X_right.empty())
        TrainDecisionNode(X_right, y_right, tree[2]);
}


// ============================
// 予測（深さ2決定木）
// ============================
// tree : 学習済みの決定木（親＋子2つ）
// x    : 1つのデータ（特徴量ベクトル）
int Predict(const vector<TreeNode>& tree, const vector<double>& x)
{
    // ===== 親ノードで分岐 =====
    if (x[tree[0].feature_id] <= tree[0].threshold) {

        // ===== 左側の枝 =====
        if (x[tree[1].feature_id] <= tree[1].threshold)
            return tree[1].left_class_id;   // 左子ノードの左側
        else
            return tree[1].right_class_id;  // 左子ノードの右側
    }
    else {

        // ===== 右側の枝 =====
        if (x[tree[2].feature_id] <= tree[2].threshold)
            return tree[2].left_class_id;   // 右子ノードの左側
        else
            return tree[2].right_class_id;  // 右子ノードの右側
    }
}


// ============================
// ラベル表示
// ============================
string LabelName(int label)
{
    if (label == 1)
        return "可溶性を示す";
    else
        return "可溶性を示さない";
}


// ============================
// 決定木表示（深さ2）
// ============================
void PrintDecisionTree(
    const vector<TreeNode>& tree,
    const vector<string>& feature_name)
{
    cout << endl;
    cout << "========== Decision Tree ==========" << endl;

    // ===== 親ノード =====
    cout << "[Root Node]" << endl;

    cout << "Feature : "
         << feature_name[tree[0].feature_id] << endl;

    cout << "Threshold : "
         << tree[0].threshold << endl;

    // ===== 左側 =====
    cout << endl;

    cout << "  ├─ NO (<= threshold)" << endl;

    cout << "  │   Feature : "
         << feature_name[tree[1].feature_id] << endl;

    cout << "  │   Threshold : "
         << tree[1].threshold << endl;

    cout << "  │   ├─ NO  : "
         << LabelName(tree[1].left_class_id) << endl;

    cout << "  │   └─ YES : "
         << LabelName(tree[1].right_class_id) << endl;

    // ===== 右側 =====
    cout << endl;

    cout << "  └─ YES (> threshold)" << endl;

    cout << "      Feature : "
         << feature_name[tree[2].feature_id] << endl;

    cout << "      Threshold : "
         << tree[2].threshold << endl;

    cout << "      ├─ NO  : "
         << LabelName(tree[2].left_class_id) << endl;

    cout << "      └─ YES : "
         << LabelName(tree[2].right_class_id) << endl;

    cout << "===================================" << endl;
}


// ============================
// 評価
// ============================
// TP: 正しく1と予測
// FP: 0なのに1と予測
// FN: 1なのに0と予測
// TN: 正しく0と予測
void Evaluation(
    const vector<TreeNode>& tree,
    const vector<vector<double>>& X,
    const vector<int>& y)
{
    int TP = 0, FP = 0, FN = 0, TN = 0; // 混同行列の各要素を初期化

    // ===== 全テストデータを評価 =====
    for (int i = 0; i < X.size(); i++) {

        int pred = Predict(tree, X[i]); //予測ラベル

        // ===== 結果を分類 =====
        if (pred == 1 && y[i] == 1) TP++;       // 正しく1
        else if (pred == 1 && y[i] == 0) FP++;  // 誤って1
        else if (pred == 0 && y[i] == 1) FN++;  // 誤って0
        else TN++;                              // 正しく0
    }

    double acc = (double)(TP + TN) / X.size(); //正解率

    // ===== Precision（適合率）=====
    double prec;

    if (TP + FP == 0)
        prec = 0;   // 分母0対策
    else
        prec = (double)TP / (TP + FP);

    // ===== Recall（再現率）=====
    double rec;

    if (TP + FN == 0)
        rec = 0;
    else
        rec = (double)TP / (TP + FN);

    // ===== F1スコア =====
    double f1;

    if (prec + rec == 0)
        f1 = 0;
    else
        f1 = 2 * prec * rec / (prec + rec);

    // ===== 結果表示 =====
    cout << "Accuracy : " << acc << endl;
    cout << "Precision: " << prec << endl;
    cout << "Recall   : " << rec << endl;
    cout << "F1-score : " << f1 << endl;

    // ===== 混同行列表示 =====
    cout << endl;

    cout << "Confusion Matrix" << endl;

    cout << "TP: " << TP << " FP: " << FP << endl;
    cout << "FN: " << FN << " TN: " << TN << endl;
}


// ============================
// main関数（全体の流れ）
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

    // ===== 学習用・テスト用データの入れ物 =====
    vector<vector<double>> training_dataset, test_dataset;

    vector<int> training_labels, test_labels;

    // ===== データ分割 =====
    DivideDataset(dataset, labels,
        training_dataset, training_labels,
        test_dataset, test_labels, 0.2);

    // ===== 深さ2の決定木（3ノード）=====
    // decision_tree[0] : 親ノード（最初の分岐）
    // decision_tree[1] : 左側の子ノード
    // decision_tree[2] : 右側の子ノード
    vector<TreeNode> decision_tree(3);

    // ===== 学習 =====
    TrainDecisionTree(training_dataset,
                      training_labels,
                      decision_tree);

    // ===== 評価 =====
    Evaluation(decision_tree,
               test_dataset,
               test_labels);

    // ===== 決定木表示 =====
    PrintDecisionTree(decision_tree, feature_name);

    return 0;
}