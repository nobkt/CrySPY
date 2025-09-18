# 分子性結晶生成例

このディレクトリには、cryspy-mol-crystallスクリプトで使用できる分子のサンプルXYZファイルが含まれています。

## サンプル分子

- `water.xyz`: 水分子 (H2O)
- `methane.xyz`: メタン分子 (CH4)
- `ammonia.xyz`: アンモニア分子 (NH3)

## 使用例

```bash
# 水分子の結晶生成
cryspy-mol-crystal water.xyz

# 複数分子の混合結晶
cryspy-mol-crystal water.xyz methane.xyz --nmol 2 2

# 3種類の分子の混合結晶
cryspy-mol-crystal water.xyz methane.xyz ammonia.xyz --nmol 1 1 2
```

## 出力例

実行すると、以下のようなファイルが生成されます：
- `optimized_crystal.cif`: 最適化された結晶構造
- `log_mol_crystal`: 実行ログ