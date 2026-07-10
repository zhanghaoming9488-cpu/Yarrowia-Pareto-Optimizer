# YaliOpt: Multi-Objective Codon Optimizer for *Yarrowia lipolytica* 🧬💻

> An advanced sequence design pipeline driven by Multi-omics Data and Non-dominated Sorting Genetic Algorithm II (NSGA-II).

![Python](https://img.shields.io/badge/Python-3.8%2B-blue.svg)
![License](https://img.shields.io/badge/License-MIT-green.svg)
![Algorithm](https://img.shields.io/badge/Algorithm-NSGA--II-orange.svg)

## 💡 Overview

**YaliOpt** is an industrial-grade bioinformatics pipeline designed to resolve the bottleneck of heterologous protein expression in *Yarrowia lipolytica*. 

Instead of traditional single-objective CAI (Codon Adaptation Index) maximization, this project utilizes a **Dual-track Strategy** (Ramp vs. Body) and a **Multi-objective Pareto Genetic Algorithm** to balance:
1. **Translation Initialization Speed** (5' Ramp CAI)
2. **Global Translation Efficiency** (Body CAI)
3. **Translational Toxicity Avoidance** (Codon Pair Score Penalty)

## 🚀 Quick Start (Demo)

Want to see the Pareto algorithm in action without downloading gigabytes of raw sequencing data? We have pre-calculated the biological rules for you!

```bash
# 1. Clone the repository
git clone [https://github.com/zhanghaoming9488-cpu/YaliOpt.git](https://github.com/zhanghaoming9488-cpu/YaliOpt.git)
cd YaliOpt

# 2. Run the Core Pareto Optimization Engine
python 03_Core_Algorithm/09_NSGA2_pareto_optimizer.py
