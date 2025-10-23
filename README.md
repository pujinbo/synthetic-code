# synthetic-code


# GitHub 仓库简介（中/英）

## 中文简介
本仓库聚焦**合成控制法（Synthetic Control Method, SCM）** 及其拓展方法的代码实现与应用指南，旨在为因果推断领域的研究者、学生及实务工作者提供一站式技术资源。  

仓库核心内容涵盖三部分：  
1. **基础实现**：包含经典SCM的核心代码（如Abadie经典框架），支持面板数据的处理、权重估计与反事实构建，附详细注释与示例数据集；  
2. **前沿拓展**：更新近年主流改进方法，重点覆盖溢出效应缓解（如迭代合成法、预设结构法）、多处理组SCM、时间异质性SCM等热点技术；  
3. **应用工具**：提供学术论文复现案例（如政策评估、区域经济研究）。  

代码以Python（Pandas、Scikit-learn）和R（Synth包、CausalImpact）为主，兼顾易用性与可扩展性，可直接适配不同研究场景的数据格式。欢迎星标（Star）关注更新，也期待通过Issues或Pull Request交流改进建议，共同推动SCM方法的实践落地。


## English Introduction
This repository focuses on the code implementation and application guidelines of the **Synthetic Control Method (SCM)** and its extensions, aiming to provide a one-stop technical resource for researchers, students, and practitioners in the field of causal inference.  

The core content of the repository includes three parts:  
1. **Basic Implementation**: Core code for classical SCM (e.g., Abadie’s framework), supporting panel data processing, weight estimation, and counterfactual construction, with detailed comments and sample datasets;  
2. **Cutting-Edge Extensions**: Updated with mainstream improved methods in recent years, focusing on spillover effect mitigation (e.g., Iterative SCM, Preset Structure Method), multi-treatment SCM, and time-heterogeneous SCM;  
3. **Application Tools**: Functions for result visualization (e.g., placebo test plots, dynamic effect plots), code for robustness checks (e.g., permutation tests, weight sensitivity analysis), and reproducible cases of academic papers (e.g., policy evaluation, regional economic research).  

Codes are mainly written in Python (Pandas, Scikit-learn) and R (Synth package, CausalImpact), balancing usability and scalability, and can be directly adapted to data formats in different research scenarios. Feel free to star the repository for updates, and we welcome suggestions via Issues or Pull Requests to jointly promote the practical application of SCM methods.

---
## 仓库代码对应文献
| File Name | Authors | Year | Title | DOI |
|-----------|------|------|----------|-----|
| hcw-data | Hsiao et al. | 2012 | A Panel Data Approach for Program Evaluation -- Measuring the Benefits of Political and Economic Integration of Hong Kong with Mainland China | 10.1002/jae.1230 |
| Abadie-2015 | Abadie et al. | 2015 | Comparative Politics and the Synthetic Control Method | 10.1111/ajps.12116 |
| gsynth-master | Xu | 2017 | Generalized Synthetic Control Method: Causal Inference with Interactive Fixed Effects Models | 10.1017/pan.2016.2 |
| robust control-main | Amjad et al. | 2017 | Robust Synthetic Control | 10.48550/arXiv.1711.06940 |
| augsynth | Ben-Michael et al. | 2018 | The Augmented Synthetic Control Method | 10.3386/w28885 |
| Firpo-2018 | Firpo and Possebom | 2018 | Synthetic Control Method: Inference, Sensitivity Analysis and Confidence Sets | 10.1515/jci-2016-0026 |
| bscm-master | Kim | 2020 | Bayesian Synthetic Control Methods | 10.1177/0022243720936230 |
| Esposti-2020 | Esposti et al. | 2020 | Can synthetic controls improve causal inference in interrupted time series evaluations of public health interventions? | 10.1093/ije/dyaa152 |
| scul-master | Hollingsworth and Wing | 2020 | Tactics for design and inference in synthetic control studies: An applied example using high-dimensional data | 10.2139/ssrn.3592088 |
| Chernozhukov-2021-JASA | Chernozhukov et al. | 2021 | An Exact and Robust Conformal Inference Method for Counterfactual and Synthetic Controls | 10.1080/01621459.2021.1920957 |
| Distributional conformal prediction | Chernozhukov et al. | 2021 | Distributional conformal prediction | 10.1073/pnas.2107794118 |
| Ferman-2021-JASA | Ferman | 2021 | On the Properties of the Synthetic Control Estimator with Many Periods and Many Controls | 10.48550/arXiv.1906.06665 |
| Ferman-2021-QE | Ferman and Pinto | 2021 | Synthetic controls with imperfect pretreatment fit | 10.3982/QE1596 |
| Fernandez 2021 | Fernandez et al. | 2021 | Low-rank approximations of nonseparable panel models | 10.1093/ectj/utab007 |
| Imai-2021 | Imai et al. | 2021 | Matching Methods for Causal Inference with Time-Series Cross-Sectional Data | 10.1111/ajps.12685 |
| MASC | Kellogg et al. | 2021 | Combining Matching and Synthetic Control to Tradeoff Biases From Extrapolation and Interpolation | 10.1080/01621459.2021.1979562 |
| pensynth | Abadie and L'Hour | 2021 | A Penalized Synthetic Control Estimator for Disaggregated Data | 10.1080/01621459.2021.1971535 |
| RPCA | Bayani | 2021 | Robust PCA Synthetic Control | 10.48550/arXiv.2108.12542 |
| scpi-main | Cattaneo et al. | 2021 | Prediction Intervals for Synthetic Control Methods | 10.1080/01621459.2021.1979561 |
| synthdid-sdid-paper | Arkhangelsky et al. | 2021 | Synthetic Difference in Differences | 10.1257/aer.20190159 |
| treebased synthetic controls-main | Mühlbach and Nielsen | 2021 | Tree-based synthetic control methods: Consequences of relocating the US embassy | 10.48550/arXiv.1909.03968 |
| synthetic learner | Viviano and Bradic | 2022 | Synthetic Learner: Model-free inference on treatments over time | 10.48550/arXiv.1904.01490 |
| CAPPMx-main | Chandra et al. | 2023 | Bayesian Nonparametric Common Atoms Regression for Generating Synthetic Controls in Clinical Trials | 10.48550/arXiv.2201.00068 |
| Cattaneo-2023 | Cattaneo et al. | 2023 | Uncertainty Quantification in Synthetic Controls with Staggered Treatment Adoption | 10.48550/arXiv.2210.05026 |
| Klinenberg-2023 | Klinenberg | 2023 | Synthetic Control with Time Varying Coefficients A State Space Approach with Bayesian Shrinkage | 10.1080/07350015.2022.2102025 |
| Replication DSC ECMA | Gunsilius | 2023 | DISTRIBUTIONAL SYNTHETIC CONTROLS | 10.3982/ECTA18260 |
| shixu-2023 | Shi et al. | 2023 | Theory for Identification and Inference with Synthetic Controls: A Proximal Causal Inference Framework | 10.48550/arXiv.2108.13935 |
| synthetic-combinations-main | Agarwal et al. | 2023 | Synthetic Combinations: A Causal Framework for Combinatorial Interventions | 10.48550/arXiv.2303.14226 |
| DR Proximal SC-main | Qiu et al. | 2024 | Doubly robust proximal synthetic controls | 10.48550/arXiv.2210.02014 |
| Liu-2024-AJPS | Liu et al. | 2024 | A Practical Guide to Counterfactual Estimators for Causal Inference with Time-Series Cross-Sectional Data | 10.1111/ajps.12723 |
| Spillover-SCM-main | Melnychuk | 2024 | Synthetic Controls with Spillover Effects: A Comparative Study | 10.48550/arXiv.2405.01645 |
| Ylmaz-2024-J0M | Yılmaz et al. | 2024 | Causal inference under selection on observables in operations management research: Matching methods and synthetic controls | 10.1002/joom.1318 |
