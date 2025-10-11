# cryspy-mol-crystal: 分子性結晶CIF生成スクリプト

分子のxyz座標から、その分子で構成される最適化分子性結晶のCIFファイルを生成するスクリプトです。

## 機能

- XYZファイルから分子データを読み込み
- pyxtalを使用して分子性結晶構造をランダム生成
- ASEを使用した構造最適化
- CIFフォーマットでの結晶構造出力
- **結晶構造情報の自動出力**（分子数、空間群番号、密度）
- 複数分子種の対応
- カスタマイズ可能なパラメータ
- 効率的な再試行メカニズムで構造生成を高速化

## インストール

CrySPYがインストールされていれば、`cryspy-mol-crystal`コマンドが利用できます。

```bash
pip install csp-cryspy
```

## 使用方法

### 基本的な使用方法

```bash
# 水分子から分子性結晶を生成
cryspy-mol-crystal water.xyz

# 出力ファイル名を指定
cryspy-mol-crystal water.xyz -o water_crystal.cif

# 複数の分子種から混合結晶を生成
cryspy-mol-crystal water.xyz methane.xyz --nmol 2 4
```

**出力ファイル:**
- `<output>.cif`: 結晶構造のCIFファイル
- `<output>_info.txt`: 結晶構造の詳細情報（分子数、空間群番号、密度など）
```

### オプション

- `-o, --output`: 出力CIFファイル名 (default: optimized_crystal.cif)
- `-n, --nstruct`: 生成する結晶構造の数 (default: 1)
- `--nmol`: 各分子の数（分子種ごとに指定、未指定時は自動設定）
- `--spgnum`: 空間群番号のリスト (default: 一般的な分子結晶の空間群)
- `--vol-factor`: 体積因子 (default: 1.1)
- `--fmax`: 構造最適化の力の閾値 (eV/Å) (default: 0.05)
- `--steps`: 構造最適化の最大ステップ数 (default: 1000)
- `--calculator`: 構造最適化に使用する計算機 (default: EMT)
- `--mindist-factor`: 最小距離の因子 (default: 0.8)
- `--no-optimization`: 構造最適化を行わない
- `--debug`: デバッグモード

### 使用例

#### 1. 水分子の結晶生成

```bash
# water.xyzファイルを用意
cat > water.xyz << EOF
3
Water molecule
O   0.000000   0.000000   0.000000
H   0.757000   0.586000   0.000000
H  -0.757000   0.586000   0.000000
EOF

# 結晶生成
cryspy-mol-crystal water.xyz --debug
```

#### 2. メタン分子の結晶生成（最適化なし）

```bash
# methane.xyzファイルを用意
cat > methane.xyz << EOF
5
Methane molecule
C   0.000000   0.000000   0.000000
H   0.627000   0.627000   0.627000
H  -0.627000  -0.627000   0.627000
H  -0.627000   0.627000  -0.627000
H   0.627000  -0.627000  -0.627000
EOF

# 結晶生成（最適化なし）
cryspy-mol-crystal methane.xyz --no-optimization -o methane_crystal.cif
```

#### 3. 水とメタンの混合結晶生成

```bash
# 混合結晶生成（水2個、メタン3個）
cryspy-mol-crystal water.xyz methane.xyz --nmol 2 3 -o mixed_crystal.cif --fmax 0.1 --steps 100
```

#### 4. 複数の結晶構造生成

```bash
# 3つの異なる結晶構造を生成
cryspy-mol-crystal water.xyz -n 3 -o water_crystals.cif
```

#### 5. 特定の空間群での生成

```bash
# P1 (空間群1) とP21/c (空間群14) のみで生成
cryspy-mol-crystal water.xyz --spgnum 1 14 -o specific_sg_crystal.cif
```

## XYZファイルフォーマット

XYZファイルは以下の形式である必要があります：

```
原子数
コメント行
原子種 x座標 y座標 z座標
原子種 x座標 y座標 z座標
...
```

例：
```
3
Water molecule
O   0.000000   0.000000   0.000000
H   0.757000   0.586000   0.000000
H  -0.757000   0.586000   0.000000
```

## 注意事項

- 構造最適化にはEMT計算機を使用します（経験的ポテンシャル、定性的な結果のみ）
- より精密な結果が必要な場合は、生成されたCIFファイルを他の量子化学計算ソフトで最適化してください
- 大きな分子や複雑な分子では結晶生成に時間がかかる場合があります
- 生成される結晶構造は初期推定であり、実験的に安定とは限りません
- 特定の空間群での構造生成が困難な場合、自動的にスキップされます（5回失敗後）
- デバッグモード（`--debug`）を使用すると、すべての試行の詳細が表示されます

## 出力ファイルの詳細

### CIFファイル (`<output>.cif`)
標準的なCIF形式で結晶構造を出力します。コメント行に密度と体積の情報も含まれます。

### 構造情報ファイル (`<output>_info.txt`)
以下の情報が人間が読みやすい形式で出力されます：
- 単位格子あたりの分子数（各分子種ごとの内訳も含む）
- 空間群番号と記号
- 密度（g/cm³）
- 体積（Å³）

例：
```
# Structure 1
Number of molecules per unit cell: 4
  Molecules per type: (4,)
Space group number: 2
Space group symbol: P-1
Density: 0.6011 g/cm³
Volume: 199.0817 Å³
============================================================
```

## 依存関係

- pymatgen
- pyxtal
- ASE (Atomic Simulation Environment)
- numpy

## ログファイル

実行時に以下のログファイルが生成されます：
- `log_mol_crystal`: 実行ログ