# `codon` — 序列分析与优化

计算每个现代 mRNA 设计模型（CodonBERT、RiboDecode、LinearDesign、mRNABERT）
都会消费的经典密码子使用特征。

## 用法

```bash
mrnavax codon --sequence mrnavax/examples/cas9.fasta
mrnavax codon --sequence mrnavax/examples/cas9.fasta --optimize
```

## 输出 schema

```json
{
  "n_codons": 210,
  "cai": 0.6954,                        // Codon Adaptation Index (0–1)
  "gc_percent": 37.3,                   // overall GC%
  "rare_codon_fraction": 0.0429,        // fraction of codons below 0.10 freq
  "cpg_obs_exp": 1.1331,                // CpG O/E ratio — proxy for innate immunity
  "most_common_codons": [["GAT", 12], ...],
  "rare_codons": ["CTA", "TTA"],
  "gc_window_stddev": 4.75               // rolling GC stddev (translation speed proxy)
}
```

## `--optimize`

运行一轮贪心式的同义密码子替换，在保持 GC% 处于 45–60% 区间的条件下，
最大化逐密码子的使用频率。

在人类细胞中表达一个细菌基因时，预期 CAI 的绝对提升约为 0.2
（例如 Cas9 示例从 0.7 提升至 0.93）。

这是任何现代密码子模型都应超越的**经典基线**。它的设计刻意保持简单，
便于你把更复杂的模型（CodonBERT、RiboDecode）接入到相同的输入/输出形状中进行对比。

## 多目标优化（知识注入式损失）

`multi-objective` 后端实现了来自 arXiv:2505.23862（RNop，*mRNA Design and Optimization with Deep Knowledge-Infused Approach*，2026 年 8 月）的 SOTA 模式：将每个生物学目标分解为各自的评分，使用者可以看到**每个**同义密码子替换是由**哪一个**分量驱动的。这正是顶级项目（RNop / LinearDesign / RiboDecode）共同汇聚的「不可能三角形」（保真度、多目标、效率）解法：保真度（仅同义替换）、多目标明确（加权损失向量）、效率以纯 Python 标准库达成。

五个归一化分量：

| 分量 | 范围 | 衡量内容 |
|---|---|---|
| `cai` | (0, 1] | Codon Adaptation Index — 转译效率 |
| `gc_score` | [0, 1] | GC% 是否落入目标区间 — 转录稳定性之代理指标 |
| `cpg_score` | [0, 1] | CpG Obs/Exp 规避 — 先天免疫之代理指标 |
| `rare_run_penalty` | [0, 1] | 最长稀有密码子连续长度 — 核糖体停滞之代理指标 |
| `structure_proxy` | [0, 1] | GC 滑动窗口标准差的一致性 — 局部结构之代理指标 |

整体分数为各分量的加权和（`rare_run_penalty` 以减项处理），夹至 `[0, 1]`。默认权重遵循 SOTA 典型模式：CAI 为主导信号（~40%），其次为结构代理（~25%），再来是 GC 与 CpG，稀有密码子连续长度计为软性惩罚。

```bash
mrnavax codon --sequence gfp.fasta \
    --optimize --backend multi-objective
```

输出包含 `before` 与 `after` 各分量的分解，加上 `overall_before`、`overall_after`、`improvement` 与 `n_changes`，使各生物轴的贡献完全可审计。

**自定义损失向量**：`MultiObjectiveConfig` 暴露全部 5 个权重（以及 GC 区间、稀有阈值、结构窗口）。在 Python 中可自行配置：

```python
from mrnavax.codon_multi_objective import (
    MultiObjectiveConfig, multi_objective_optimize,
)

cfg = MultiObjectiveConfig(
    cai_weight=0.30,        # 略降 CAI 比重，提高结构权重
    gc_weight=0.20,
    cpg_weight=0.15,
    rare_run_weight=0.25,   # 强调避免稀有密码子连续
    structure_weight=0.10,
    target_gc_min=50.0,     # 更紧的 GC 区间
    target_gc_max=58.0,
)
result = multi_objective_optimize(cds, cfg)
print(result.improvement, result.n_changes)
```

## RNop 参考文献

Gong, Z., Jiang, Z., Gao, W., Wang, Y., Cai, Z., Deng, Z., Ma, L.
(2026). *mRNA Design and Optimization with Deep Knowledge-Infused
Approach.* arXiv:2505.23862v2。「不可能三角形」的论述（保真度、多目标、效率）与知识注入式损失分解为 `multi-objective` 后端分量分解的概念基础。
