# `utr-design` — 耦合 5'UTR + CDS + 3'UTR 以最大化表达量

mrnavax 工具组中的**第 10 个工具**。将 v0.27.0 的 UTR 上下文评分器
（Kozak + 3'UTR 质量）与 v0.25.0 的多目标 CDS 优化器
（RNop 知识注入式损失模式）耦合，针对给定的蛋白质氨基酸序列
挑选最佳组合。

## 为何称为「耦合」？

v0.26.0 的 `construct` 工具从固定的 UTR 模板组出 5'UTR + CDS + 3'UTR + poly-A
构建——对组装有帮助，但 UTR 的选择是写死的。v0.27.0 的 UTR 上下文评分器
可按表达量影响为 UTR 排序，但无法帮忙挑选。

`utr-design` 就是**桥梁**：在 12 个候选 UTR 组合（4 个 5'UTR 变体 × 3 个 3'UTR 变体）
上进行有界网格搜索，逐一用 v0.27.0 UTR 上下文评分器评分，
CDS 用 v0.25.0 多目标优化器优化一次，再以最高**联合分数**
（0.4 × UTR 上下文 + 0.6 × CDS）选出组合。

## 快速上手

```python
from mrnavax.utr_designer import UTRDesignConfig, design_utr_aware_cds

cfg = UTRDesignConfig(cds="MVSKGEELFTGV")  # 12-AA eGFP 片段
result = design_utr_aware_cds(cfg)

print(f"Chosen UTR5: {result.utr5}")
print(f"Chosen UTR3: {result.utr3}")
print(f"Optimized CDS: {result.cds_dna}")
print(f"Combined score: {result.combined_score:.3f}")
```

## CLI

```bash
mrnavax utr-design --protein "MVSKGEELFTGV"
mrnavax utr-design --protein my_protein.fasta --prefer-kozak 0.85
mrnavax utr-design --protein my_protein.fasta --library minimal
```

## 候选 UTR 库

默认 `library="all"` 评估 12 种组合：

|| 5'UTR 变体 | 来源 | Kozak 分数 |
||---|---|---|
|| `strong_kozak` | 9-nt 片段 `GCCACCAAT`（R=A 嘌呤） | ~0.9 |
|| `moderate_kozak` | 9-nt 片段 `GCATGG`（变体） | ~0.5 |
|| `weak_kozak` | 9-nt 片段 `AAAAAA..`（无 GCC） | ~0.1 |
|| `mrna1273_like` | mRNA-1273 已发表的 5'UTR | ~0.5 |

|| 3'UTR 变体 | 来源 | ARE 负担 |
||---|---|---|
|| `short_constitutive` | 32 nt GC 平衡，无 ARE | ~1.0 |
|| `mrna1273_like` | mRNA-1273 已发表的 3'UTR | ~0.7 |
|| `long_conservative` | 60 nt GC 平衡 | ~1.0 |

`library="minimal"` 设置仅使用 strong-Kozak + short-constitutive 对（1 种组合）。

## 自定义阈值

默认情况下，设计师要求 Kozak ≥ 0.7 且 3'UTR ≥ 0.5。若没有任何组合同时
满足两个阈值，则回退选取整体排名最高的组合。

```python
# 对 Kozak 严格、对 3'UTR 宽松
cfg = UTRDesignConfig(
    cds="MVSKGEELFTGV",
    prefer_kozak=0.85,
    prefer_utr3=0.3,
)
```

## 选择 CDS 后端

默认情况下，CDS 以多目标后端（v0.25.0 RNop 知识注入式损失）进行优化。
`basic` 后端亦可使用：

```python
cfg = UTRDesignConfig(cds="MVSKGEELFTGV", backend="basic")
```

推荐多目标后端，是因为它在单一加权损失中同时考虑
CAI + GC + CpG + 稀有连读 + 结构代理，给出比 CAI+GC 更细致的 CDS 选择。

## 输出

结果是一个 `UTRDesignResult` 数据类，内含组装最终构件所需的全部信息：

```python
@dataclass
class UTRDesignResult:
    utr5: str
    utr3: str
    cds_dna: str
    protein: str
    utr_context: UTRContextResult
    cds_score: float
    combined_score: float
    candidates_evaluated: int
    ranking: list[tuple[float, str, str]]  # (combined_score, utr5_name, utr3_name)
```

使用 `result.utr5 + result.cds_dna + result.utr3` 作为完整 mRNA 构建中
「CDS 两侧接 UTR」的部分。

## 未优化的部分

- **5'UTR / 3'UTR 的序列内容**：我们从手工策展的库（4 + 3 = 7 个变体）
  中采样。若要做序列层级的 UTR 优化（这是另一个问题——在可微分
  表达量模型上做梯度式设计），请使用专用工具如 UTR-Function 或
  5UTR-Tailor。库方法对不确定性更诚实，也能避免过拟合。
- **跨 UTR 配对的 CDS 逐位置选择**：CDS 只优化一次，之后在 12 种 UTR 配对
  中重复使用。这在 CDS 与 UTR 选择相互独立的意义下是正确的，但也意味着
  我们并未对 (CDS, UTR) 做联合搜索。实际上主导信号来自 CDS 本身，
  所以这样做没有问题。
- **Poly-A 尾长度**：完整的 poly-A 组装请见 `construct` 工具。
  `utr-design` 专注于 UTR × CDS 的耦合。