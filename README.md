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
