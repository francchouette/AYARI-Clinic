# Source data and scientific provenance

本GLBは、Janelia Research Campus の OpenOrganelle / COSEM 公開データセット
`jrc_hela-2` を形態参照として使用したハイブリッドモデルです。

## Primary measured source

- Dataset: **Automatic segmentation of organelles in isotropic electron microscopy data of wild-type interphase HeLa cell (`jrc_hela-2`)**
- Sample: wild-type interphase HeLa cell, ATCC CCL-2
- Modality: FIB-SEM
- Native voxel size: `4 × 4 × 5.24 nm`
- Dataset page: <https://openorganelle.janelia.org/datasets/jrc_hela-2>
- DOI and downloadable segmentations: <https://doi.org/10.25378/janelia.13108343>
- License: **CC BY 4.0**
- Attribution: Janelia Research Campus / OpenOrganelle COSEM, dataset `jrc_hela-2`

GLB生成では、Web表示用のポリゴン予算に合わせて、公開N5の細胞foreground
`masks/foreground/s3` とplasma-membrane label `pm_seg/s4`、および細胞小器官の
`s4` labelsから表面を抽出しています。これらは同じ `64 × 64 × 83.84 nm`
sampling gridです。元データの出自を隠さないため、GLBのasset extrasと各部位
node extrasにもdataset、DOI、source label、scale、source instance ID、sampling
解像度を格納しています。

| GLB node | N5 source label | 扱い |
|---|---|---|
| `GEO_cell` | `masks/foreground` (`s3`), `pm_seg` instance 2 (`s4`) | 実測対象細胞の境界をplasma-membrane labelで補強した表面 |
| `GEO_nucleus` | `nucleus_seg`, `chrom_seg` | 実測セグメンテーション表面 |
| `GEO_mitochondria_01…07` | `mito_seg`, `mito-mem_seg` | 実測個体外形＋同一IDの膜／クリステ |
| `GEO_chromosome_*`, `GEO_telomere_*` | なし | 教育用オーバーレイ（非実測） |
| `GEO_chromatin_fiber`, `GEO_nucleosome_*`, `GEO_dna` | なし | 核内DNAのパッケージング階層を拡大した教育用オーバーレイ（非実測） |

公開foreground maskには中央の対象細胞と、基底面で接触する隣接細胞の一部が
含まれます。対象細胞は最大の中央核をseedとしてY slice方向にforeground成分を
追跡し、成分が隣接細胞や取得境界へ接続するsliceだけseeded watershedで接触面を
分離しています。したがって外周は実測foreground境界、接触分離面は再構成境界です。
この選択方法は `GEO_cell.extras.targetCellSelection` にも記録しています。

## Display-space normalization

元のHeLa細胞は培養基板上で扁平な実測形状です。本アセットでは、依頼された
直径約2 mのカットアウェイ構図と部位選択UIに合わせ、実測細胞膜を
`2.00 × 1.84 × 1.72 m` の表示空間へ非一様正規化しています。核は実測
トポロジーを保ったまま指定の約0.60 m幅へ非一様正規化し、ミトコンドリアは
各個体をPCA整列後に一様スケールして構図上へ再配置しています。この変換は
形態比較・距離計測・体積計測には使用できません。

`jrc_hela-2` は**間期**細胞です。完全な核の内部に凝縮したX型染色体が存在する
状態を同一時点の実測像として扱うのは正確ではないため、X型染色体、末端の
テロメアキャップ、DNAパッケージング階層は明示的に教育用オーバーレイとして
います。DNAは核前面で `GEO_chromosome_01` へ接続し、染色体→簡略30 nm相当
クロマチン繊維→ヒストン八量体と約1.65回巻きDNAからなるヌクレオソーム→
右巻きDNA二重らせん、という拡大関係を構造で示します。実寸2 nm、11 nm、
30 nmはnode extrasに保存していますが、各形状は視認性のため同一縮尺ではありません。

## Reproduction

```bash
python tools/fetch_openorganelle_cell.py
python tools/build_cell_glb.py
python tools/validate_cell_glb.py \
  assets/cell/ayari_cell_cutaway.glb \
  --report assets/cell/model_report.json
```

取得したN5キャッシュは `.cache/` 配下に保存され、Gitへは追加されません。
再生成スクリプトが公開S3から必要なlabel blockを取得します。

## Use limitation

本アセットは研究データを参照した教育・可視化用モデルです。診断、病理判定、
定量的な細胞計測、治療計画、手術計画には使用できません。
