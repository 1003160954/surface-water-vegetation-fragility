# Surface water permanence and vegetation fragility

Project repository: https://github.com/1003160954/surface-water-vegetation-fragility

可复核的水体背景与植被干热风险分析 / Reproducible analysis of aquatic context and vegetation risk amplification across China, 2018–2023.

本仓库公开**代码、数据来源说明、派生结果表和可视化**。论文与审稿回复工作稿不在公开范围内；大规模原始导出表不直接存入Git。结果来自真实分析，图形由数值表与绘图程序生成。使用过的AI辅助另行如实披露；“可复现”不等于“没有使用AI”，也不等于没有研究局限。

## 先看研究结果

问题是：夏天又热又干时，永久水体背景附近的植被，风险放大是否较弱？“风险放大”是高干热组事件比例减低干热组事件比例。“事件”是植被指数相对自身近期季节水平偏低，并非死亡或永久植被损失。

- 主样本236,995条网格年、48,843个0.1°网格。纠正后的永久、季节性、其他背景分别有2,549、15,262、31,032个网格。
- 永久背景的EVI、LAI、组合事件风险放大分别为0.1292、0.1408、0.1600；其他背景分别为0.1901、0.1972、0.2209。原始事件计数与除法运算见[整数计数表](reviewer_evidence/risk_numerators_and_arithmetic.csv)。
- 对其他背景的三个差异，在已测协变量调整和排除主要农业覆盖后均保留网格聚类区间支持；永久相对季节性水体的证据以EVI最清楚，LAI不确定，组合事件依目标样本和不确定性设置而变化。
- 1°和2°空间块下，三个其他背景对比和EVI永久性对比仍有支持；5°下全部未调整对比区间包含零。全部尺度均报告，不能只挑有利结果，也不能把一个跨零区间当成“证明无效”。
- 这是**缓冲关联**，不是已经识别的因果水文作用；1km中心分类也不是格内植物的保护半径。

### 纠正分类后的完整调整结果

AFB = fragility(comparison) − fragility(permanent)。正值表示永久背景的风险放大较小；单位是绝对概率差。下表使用同一协变量标准化人群，区间为网格聚类sandwich/delta95%区间，不能冒充空间块区间。

| Comparison | Response | Adjusted AFB | 95% interval |
| --- | --- | --- | --- |
| Seasonal − permanent | EVI | 0.0378 | [0.0126, 0.0630] |
| Other − permanent | EVI | 0.0664 | [0.0421, 0.0908] |
| Seasonal − permanent | LAI | 0.0079 | [-0.0155, 0.0314] |
| Other − permanent | LAI | 0.0412 | [0.0186, 0.0638] |
| Seasonal − permanent | combined | 0.0323 | [0.0054, 0.0591] |
| Other − permanent | combined | 0.0667 | [0.0407, 0.0926] |

## 如何从数据走到结果

请按下面顺序阅读，而不必先理解全部代码。每步的输入、处理、输出、命令与限制见[逐步计算说明](docs/COMPUTATION_STEPS.md)。字段与缺失数量见[完整字段字典](reviewer_evidence/all_panel_fields_dictionary.csv)。

| step | operation | evidence | code |
| --- | --- | --- | --- |
| 1 | Original data products | docs/data_sources.csv; reviewer_evidence/all_panel_fields_dictionary.csv | gee/; docs/original_sources/upstream_seasonal_source_only.ipynb |
| 2 | Supplied raw analysis export | logs/E00_manifest.json; reviewer_evidence/raw_input_manifest.json | src/audit_pipeline.py |
| 3 | Validity and seasonal anomaly recipe | experiments/E00_original_reproduction/sample_flow.csv | docs/original_sources/upstream_seasonal_source_only.ipynb |
| 4 | Stress and event assignments | reviewer_evidence/worked_grid_year_examples.csv | src/audit_pipeline.py |
| 5 | Corrected aquatic contexts | experiments/E01_spatial_support/national_corrected_1km_contexts.csv | gee/corrected_near_contexts.py |
| 6 | Event risk and fragility | reviewer_evidence/risk_numerators_and_arithmetic.csv; reviewer_evidence/AFB_arithmetic.csv | src/audit_pipeline.py; src/analyze_corrected_1km.py |
| 7 | Measured-confounder adjustment | experiments/E01_spatial_support/corrected_1km/harmonized/ | src/harmonized_confounding.py; src/analyze_corrected_1km.py |
| 8 | Landscape and agricultural sensitivity | reviewer_evidence/corrected_covariate_balance.csv; experiments/E02_confounding/corrected_noncrop/ | src/common_landscape.py; src/reviewer_noncrop.py |
| 9 | Robustness and spatial uncertainty | experiments/E01_spatial_support/corrected_1km/robustness/ | src/robustness.py; src/analyze_corrected_1km.py |
| 10 | Quantitative visualizations | manuscript_revision/figure_source_map.csv; reviewer_evidence/figure_source_map.csv | src/revision_figures.py; src/build_reviewer_evidence.py |

### 用一条计算理解AFB

永久背景组合事件：低HDSI风险0.214610，高HDSI风险0.374590，风险放大0.159979。其他背景：0.192537、0.413478，风险放大0.220940。因此其他减永久AFB = 0.220940 − 0.159979 = 0.060961。

完整精度与整数分子/分母在[risk_numerators_and_arithmetic.csv](reviewer_evidence/risk_numerators_and_arithmetic.csv)和[AFB_arithmetic.csv](reviewer_evidence/AFB_arithmetic.csv)。不要用已四舍五入的图中数字验证末位精度。27条[网格年示例](reviewer_evidence/worked_grid_year_examples.csv)来自按固定随机种子选取的6个网格，仅说明计算，不是推断样本。

## 每类数据从哪里来

| 数据 | 本次实际用途 | 可核对来源与状态 |
| --- | --- | --- |
| JRC Global Surface Water v1.4 | 水体掩膜、中心分类、原生距离核查 | [官方目录](https://developers.google.com/earth-engine/datasets/catalog/JRC_GSW1_4_GlobalSurfaceWater)，JRC/GSW1_4/GlobalSurfaceWater；30m；长期产品1984–2021，非逐年水体可用性观测 |
| MODIS MCD12Q1 v6.1，2018 | 调整与分层使用的分类土地覆盖 | [官方目录](https://developers.google.com/earth-engine/datasets/catalog/MODIS_061_MCD12Q1)，MODIS/061/MCD12Q1；网格空间众数，不能平均类别码 |
| WorldClim v1 monthly，1960–1990 | 年平均温度与年降水背景协变量 | [官方目录](https://developers.google.com/earth-engine/datasets/catalog/WORLDCLIM_V1_MONTHLY)，WORLDCLIM/V1/MONTHLY；不是独立植被异常基线 |
| ESA WorldCover v200，2021 | 24网格内部植被距离试验的植被掩膜 | [官方目录](https://developers.google.com/earth-engine/datasets/catalog/ESA_WorldCover_v200)，ESA/WorldCover/v200；不用于全国主响应定义 |
| 提供的植被—气候分析导出表 | EVI/LAI、夏季气候指标、有效性标记、原背景变量 | [输入身份与SHA256](reviewer_evidence/raw_input_manifest.json)；它是分析导出表，不是原始卫星影像。原MODIS/ERA5提取部分精确collection版本与QA信息尚未全部核实，不猜填 |
| 已有climatology导出表 | 审计近期季节参考样本数 | 最多6年；[E03状态](experiments/E03_independent_baseline/status.json)。不能当成2001–2017独立参考 |

更详细的来源见[data_sources.csv](docs/data_sources.csv)和[数据访问与来源边界](docs/DATA_ACCESS.md)。水体地图请注明Source: EC JRC/Google，并遵守各产品的归属和使用条款。来源链接不是作者加工代码已完整恢复的证明。

## 怎样减少其他因素的影响

我们不以“无法避免”代替处理。已实施：本地季节标准化；地形/土地覆盖/背景气候/年份/地理项调整；有共同支持的景观分层与共同权重；排除主要耕地/镶嵌类别；网格及多尺度空间不确定性检查。[完整对照表](reviewer_evidence/confounding_controls_and_limits.csv)逐项说明能解决与不能解决什么。

在限制并共同加权的景观样本中，海拔、坡度、温度、降水的绝对描述性SMD分别约0.037、0.002、0.028、0.033。限制前后目标人群不同，不能把改善全部归因于权重，更不能把这些数字称为所有混杂都已消除。

排除2018众数土地覆盖的耕地12和耕地/自然植被镶嵌14类后，剩172,463行、36,238网格；调整模型拟合82,743个低/高HDSI记录、30,212网格。季节减永久EVI AFB为0.0480[0.0170,0.0790]，组合事件0.0410[0.0080,0.0739]，LAI0.0218[−0.0071,0.0507]。三个其他背景对比仍为正区间。见[所有排除结果](experiments/E02_confounding/corrected_noncrop/M2_AFB.csv)。这不是“无灌溉实验”；地下水、土壤、类外灌溉和动态管理仍未直接控制。

## 图表导航

### 当前纠正结果的主图

#### Fig1_workflow — 步骤与计算框架

![步骤与计算框架](revised_figures/figures/Fig1_workflow.png)

[PDF](revised_figures/figures/Fig1_workflow.pdf) · [主图来源表](revised_figures/figure_source_map.csv)

#### Fig2_JRC_attributes — JRC发生频率、季节性和离散转移类别

![JRC发生频率、季节性和离散转移类别](revised_figures/figures/Fig2_JRC_attributes.png)

[PDF](revised_figures/figures/Fig2_JRC_attributes.pdf) · [主图来源表](revised_figures/figure_source_map.csv)

#### Fig3_corrected_contexts — 全国纠正后的背景分类及组合事件空间分布

![全国纠正后的背景分类及组合事件空间分布](revised_figures/figures/Fig3_corrected_contexts.png)

[PDF](revised_figures/figures/Fig3_corrected_contexts.pdf) · [主图来源表](revised_figures/figure_source_map.csv)

#### Fig4_exposure_risks — 干热暴露的空间分布与事件风险

![干热暴露的空间分布与事件风险](revised_figures/figures/Fig4_exposure_risks.png)

[PDF](revised_figures/figures/Fig4_exposure_risks.pdf) · [主图来源表](revised_figures/figure_source_map.csv)

#### Fig5_context_fragility — 三类背景下的植被风险放大

![三类背景下的植被风险放大](revised_figures/figures/Fig5_context_fragility.png)

[PDF](revised_figures/figures/Fig5_context_fragility.pdf) · [主图来源表](revised_figures/figure_source_map.csv)

#### Fig6_corrected_contrasts — 完整未调整与调整AFB及区间

![完整未调整与调整AFB及区间](revised_figures/figures/Fig6_corrected_contrasts.png)

[PDF](revised_figures/figures/Fig6_corrected_contrasts.pdf) · [主图来源表](revised_figures/figure_source_map.csv)

#### Fig7_landscape_contrasts — 可比景观下的永久—季节性比较

![可比景观下的永久—季节性比较](revised_figures/figures/Fig7_landscape_contrasts.png)

[PDF](revised_figures/figures/Fig7_landscape_contrasts.pdf) · [主图来源表](revised_figures/figure_source_map.csv)

#### Fig8_spatial_uncertainty — 全部1°、2°、5°空间块结果

![全部1°、2°、5°空间块结果](revised_figures/figures/Fig8_spatial_uncertainty.png)

[PDF](revised_figures/figures/Fig8_spatial_uncertainty.pdf) · [主图来源表](revised_figures/figure_source_map.csv)

### 新增的纠正结果与诊断图

#### SFig12_noncrop_sensitivity

![SFig12_noncrop_sensitivity](reviewer_evidence/figures/SFig12_noncrop_sensitivity.png)

[PDF](reviewer_evidence/figures/SFig12_noncrop_sensitivity.pdf) · [数值来源](reviewer_evidence/figure_source_map.csv)

#### SFig13_corrected_balance

![SFig13_corrected_balance](reviewer_evidence/figures/SFig13_corrected_balance.png)

[PDF](reviewer_evidence/figures/SFig13_corrected_balance.pdf) · [数值来源](reviewer_evidence/figure_source_map.csv)

#### SFig14_corrected_HDSI

![SFig14_corrected_HDSI](reviewer_evidence/figures/SFig14_corrected_HDSI.png)

[PDF](reviewer_evidence/figures/SFig14_corrected_HDSI.pdf) · [数值来源](reviewer_evidence/figure_source_map.csv)

#### SFig15_corrected_thresholds

![SFig15_corrected_thresholds](reviewer_evidence/figures/SFig15_corrected_thresholds.png)

[PDF](reviewer_evidence/figures/SFig15_corrected_thresholds.pdf) · [数值来源](reviewer_evidence/figure_source_map.csv)

#### SFig16_corrected_year_omissions

![SFig16_corrected_year_omissions](reviewer_evidence/figures/SFig16_corrected_year_omissions.png)

[PDF](reviewer_evidence/figures/SFig16_corrected_year_omissions.pdf) · [数值来源](reviewer_evidence/figure_source_map.csv)

#### SFig17_sample_flow

![SFig17_sample_flow](reviewer_evidence/figures/SFig17_sample_flow.png)

[PDF](reviewer_evidence/figures/SFig17_sample_flow.pdf) · [数值来源](reviewer_evidence/figure_source_map.csv)

### 原结果复现与历史距离审计图

以下包含旧分组、原物理距离或原算法的审计结果。它们用于解释原结果怎样计算出来，**不能替代上面的纠正结果，也不能据其公里曲线推断生态作用范围**。其中SFig9–11已经是纠正1km的补充展示；具体范围见来源表。

<details>
<summary>展开全部原审计图及空间图</summary>

**Fig1_verified_workflow**

![Fig1_verified_workflow](figures/Fig1_verified_workflow.png)

[PDF](figures/Fig1_verified_workflow.pdf)

**Fig2_JRC_metrics_discrete_transition**

![Fig2_JRC_metrics_discrete_transition](figures/Fig2_JRC_metrics_discrete_transition.png)

[PDF](figures/Fig2_JRC_metrics_discrete_transition.pdf)

**Fig3_context_and_combined_events**

![Fig3_context_and_combined_events](figures/Fig3_context_and_combined_events.png)

[PDF](figures/Fig3_context_and_combined_events.pdf)

**Fig4_HDSI_and_event_risks**

![Fig4_HDSI_and_event_risks](figures/Fig4_HDSI_and_event_risks.png)

[PDF](figures/Fig4_HDSI_and_event_risks.pdf)

**Fig5_numeric_distance_bands**

![Fig5_numeric_distance_bands](figures/Fig5_numeric_distance_bands.png)

[PDF](figures/Fig5_numeric_distance_bands.pdf)

**Fig6_context_risk_and_fragility**

![Fig6_context_risk_and_fragility](figures/Fig6_context_risk_and_fragility.png)

[PDF](figures/Fig6_context_risk_and_fragility.pdf)

**Fig7_exploratory_diagnostics**

![Fig7_exploratory_diagnostics](figures/Fig7_exploratory_diagnostics.png)

[PDF](figures/Fig7_exploratory_diagnostics.pdf)

**Fig8_all_regions_with_support**

![Fig8_all_regions_with_support](figures/Fig8_all_regions_with_support.png)

[PDF](figures/Fig8_all_regions_with_support.pdf)

**SFig10_corrected_1km_spatial_blocks**

![SFig10_corrected_1km_spatial_blocks](figures/SFig10_corrected_1km_spatial_blocks.png)

[PDF](figures/SFig10_corrected_1km_spatial_blocks.pdf)

**SFig11_corrected_1km_context_map**

![SFig11_corrected_1km_context_map](figures/SFig11_corrected_1km_context_map.png)

[PDF](figures/SFig11_corrected_1km_context_map.pdf)

**SFig1_same_sample_adjustment**

![SFig1_same_sample_adjustment](figures/SFig1_same_sample_adjustment.png)

[PDF](figures/SFig1_same_sample_adjustment.pdf)

**SFig2_spatial_block_uncertainty**

![SFig2_spatial_block_uncertainty](figures/SFig2_spatial_block_uncertainty.png)

[PDF](figures/SFig2_spatial_block_uncertainty.pdf)

**SFig3_HDSI_variants**

![SFig3_HDSI_variants](figures/SFig3_HDSI_variants.png)

[PDF](figures/SFig3_HDSI_variants.pdf)

**SFig4_event_thresholds**

![SFig4_event_thresholds](figures/SFig4_event_thresholds.png)

[PDF](figures/SFig4_event_thresholds.pdf)

**SFig5_leave_analysis_year_out**

![SFig5_leave_analysis_year_out](figures/SFig5_leave_analysis_year_out.png)

[PDF](figures/SFig5_leave_analysis_year_out.pdf)

**SFig6_centroid_distance_spline**

![SFig6_centroid_distance_spline](figures/SFig6_centroid_distance_spline.png)

[PDF](figures/SFig6_centroid_distance_spline.pdf)

**SFig7_within_cell_pilot**

![SFig7_within_cell_pilot](figures/SFig7_within_cell_pilot.png)

[PDF](figures/SFig7_within_cell_pilot.pdf)

**SFig8_harmonized_and_landscape**

![SFig8_harmonized_and_landscape](figures/SFig8_harmonized_and_landscape.png)

[PDF](figures/SFig8_harmonized_and_landscape.pdf)

**SFig9_corrected_1km_permanence**

![SFig9_corrected_1km_permanence](figures/SFig9_corrected_1km_permanence.png)

[PDF](figures/SFig9_corrected_1km_permanence.pdf)

</details>

独立原生像元距离核查另见[Supplementary distance validation](revised_figures/figures/SFig18_distance_validation.png)。全部可用图、数值源、范围和脚本收录在[FIGURE_GALLERY.md](docs/FIGURE_GALLERY.md)；没有未完成实验的虚构图。

## 在本地复现

Python3.12；本项目记录的实际依赖版本见requirements.txt。大数据不在Git，设置AQUATIC_PANEL、AQUATIC_CLIMATOLOGY到作者提供的原始分析导出表。访问原始表的权限与访问本仓库的权限不同。仅有公开派生表时，可以核对算术和重画部分结果，不能声称完成从原始影像到全部模型的独立重建。

```text
python -m pip install -r requirements.txt
# Set AQUATIC_PANEL and optionally AQUATIC_CLIMATOLOGY in your shell.
python src/original_exact.py
python src/analyze_corrected_1km.py
python src/reviewer_noncrop.py
python src/build_reviewer_evidence.py
python src/revision_figures.py
python src/plot_corrected_fragility.py
python tests/validate_evidence.py
python tests/validate_public_evidence.py
```

`revision_figures.py`重建旧主图布局，其中旧Fig5距离审计仍在；`plot_corrected_fragility.py`另建当前Fig5背景比较。两种图都保留，来源不混用。完整历史复现入口与每步预期输出见COMPUTATION_STEPS.md。复现并不要求重新提交GEE导出。

## 已完成与未完成

| 项目 | 状态 |
| --- | --- |
| 原主要风险和6个AFB的精确复现 | 已完成；原500次行级bootstrap与新300次网格bootstrap区分报告 |
| 全国纠正1km中心分类 | 已完成48,843网格，17,056网格完整距离对照零分类不一致 |
| 调整、可比景观、非主要农业敏感性 | 已完成，残余混杂保留 |
| 10种HDSI、4阈值、6分析年剔除、1/2/5度块 | 已完成，全部数值与区间公开 |
| 全国完整距离幅度、全国网格内部指标 | 部分输出保留，重型任务按作者要求暂停；不得自动恢复 |
| 位移/分辨率变化后重新聚合响应和气候 | 未完成；度量掩膜移动试验不等于此项 |
| 独立2001–2017异常参考 | 未完成；软件测试不冒充真实生态实验 |
| 原植被/热产品的所有collection与QA信息 | 未全部核实；来源缺口明确标注 |
| 纠正分类后的全国机制与区域作用 | 未建立；旧输出仅作审计 |

夏季152–243与夏末213–273按16天步中点DOY划分，有213–243重叠；不能写成完全先行暴露。跨零区间不证明没有效应；正点估计不证明稳健或因果。5°空间块不确定性必须与较小尺度结果共同报告。

## 版本与文件范围

公开Git历史从筛选后的研究材料开始，不包含私有论文/回复草稿、Word文件、原始大表、凭据或旧本地仓库历史。以前未找到笔记本的初次运行记录保留时间含义；现已找回原研究和上游季节笔记本，其无历史输出的源代码归档见docs/original_sources。运行状态与文件清单见[REPRODUCIBILITY.md](docs/REPRODUCIBILITY.md)及public_manifest.csv。
