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

GLB生成では、Web表示用のポリゴン予算に合わせて、公開N5の細胞小器官 `s4`
labelsから核・クロマチン・ミトコンドリア表面を抽出しています。細胞foreground
`masks/foreground/s3` とplasma-membrane label `pm_seg/s4` は、表示用の丸みを
優先した有機エンベロープの形態参照とprovenanceに使用しています。これらは同じ
`64 × 64 × 83.84 nm` sampling gridです。元データの出自を隠さないため、GLBのasset extrasと各部位
node extrasにもdataset、DOI、source label、scale、source instance ID、sampling
解像度を格納しています。

| GLB node | N5 source label | 扱い |
|---|---|---|
| `GEO_cell` | `masks/foreground` (`s3`), `pm_seg` instance 2 (`s4`) | 実測ラベルを形態参照とし、培養時の扁平化・接触面を除いた丸みのある有機的な表示エンベロープ（非実測形状） |
| `GEO_nucleus` | `nucleus_seg`, `chrom_seg` | 実測セグメンテーション表面 |
| `GEO_mitochondria_01…07` | `mito_seg`, `mito-mem_seg` | 実測個体外形＋同一IDの膜／クリステ |
| `GEO_chromosome_*`, `GEO_telomere_*` | なし | 教育用オーバーレイ（非実測） |
| `GEO_chromatin_fiber`, `GEO_nucleosome_*`, `GEO_dna` | なし | 核内DNAのパッケージング階層を拡大した教育用オーバーレイ（非実測） |

公開foreground maskには中央の対象細胞と、基底面で接触する隣接細胞の一部が
含まれます。対象細胞は最大の中央核をseedとしてY slice方向にforeground成分を
追跡し、成分が隣接細胞や取得境界へ接続するsliceだけseeded watershedで接触面を
分離しています。この抽出結果自体は表示メッシュへ直接転写せず、sourceの選択方法と
形態参照として使用します。選択方法は `GEO_cell.extras.targetCellSelection`、表示形状の
扱いは `measured: false`、`sourceMeasuredReference: true` として記録しています。

## Display-space normalization

元のHeLa細胞は培養基板上で扁平な実測形状です。本アセットでは、直径約2 mの
カットアウェイ構図と部位選択UIに合わせ、細胞膜を低周波の非対称な膨らみを持つ
滑らかな有機形状として再構成しています。完全な球ではありませんが、培養基板による
扁平さや隣接細胞との接触分離面は表示シルエットへ持ち込みません。核は実測
トポロジーを保ったまま指定の約0.60 m幅へ非一様正規化し、ミトコンドリアは
各個体をPCA整列後に一様スケールして構図上へ再配置しています。この変換は
形態比較・距離計測・体積計測には使用できません。

`jrc_hela-2` は**間期**細胞です。完全な核の内部に凝縮したX型染色体が存在する
状態を同一時点の実測像として扱うのは正確ではないため、X型染色体、末端の
テロメアキャップ、DNAパッケージング階層は明示的に教育用オーバーレイとして
います。DNAは核内で `GEO_chromosome_01` へ接続し、染色体→簡略30 nm相当
クロマチン繊維→ヒストン八量体と約1.65回巻きDNAからなるヌクレオソーム→
右巻きDNA二重らせん、という拡大関係を構造で示します。実寸2 nm、11 nm、
30 nmはnode extrasに保存していますが、各形状は視認性のため同一縮尺ではありません。
DNA二重らせん、クロマチン繊維、3個のヌクレオソームは核の凸包内に完全包含することを
生成時とvalidatorで確認します。

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
