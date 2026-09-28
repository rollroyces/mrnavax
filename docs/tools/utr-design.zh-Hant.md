# `utr-design` — 耦合 5'UTR + CDS + 3'UTR 以最大化表現量

mrnavax 工具組中的**第 10 個工具**。將 v0.27.0 的 UTR 情境評分器
（Kozak + 3'UTR 品質）與 v0.25.0 的多目標 CDS 優化器
（RNop 知識注入式損失模式）耦合，針對給定的蛋白質胺基酸序列
挑選最佳組合。

## 為何稱為「耦合」？

v0.26.0 的 `construct` 工具從固定的 UTR 模板組出 5'UTR + CDS + 3'UTR + poly-A
構件——對組裝有幫助，但 UTR 的選擇是寫死的。v0.27.0 的 UTR 情境評分器
可依表現量影響為 UTR 排序，但無法幫忙挑選。

`utr-design` 就是**橋樑**：在 12 個候選 UTR 組合（4 個 5'UTR 變體 × 3 個 3'UTR 變體）
上做有界格點搜尋，逐一以 v0.27.0 UTR 情境評分器評分，
CDS 以 v0.25.0 多目標優化器優化一次，再以最高**聯合分數**
（0.4 × UTR 情境 + 0.6 × CDS）選出組合。

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

## 候選 UTR 函式庫

預設 `library="all"` 評估 12 種組合：

|| 5'UTR 變體 | 來源 | Kozak 分數 |
||---|---|---|
|| `strong_kozak` | 9-nt 片段 `GCCACCAAT`（R=A 嘌呤） | ~0.9 |
|| `moderate_kozak` | 9-nt 片段 `GCATGG`（變異型） | ~0.5 |
|| `weak_kozak` | 9-nt 片段 `AAAAAA..`（無 GCC） | ~0.1 |
|| `mrna1273_like` | mRNA-1273 已發表的 5'UTR | ~0.5 |

|| 3'UTR 變體 | 來源 | ARE 負擔 |
||---|---|---|
|| `short_constitutive` | 32 nt GC 平衡，無 ARE | ~1.0 |
|| `mrna1273_like` | mRNA-1273 已發表的 3'UTR | ~0.7 |
|| `long_conservative` | 60 nt GC 平衡 | ~1.0 |

`library="minimal"` 設定僅使用 strong-Kozak + short-constitutive 對（1 種組合）。

## 自訂閾值

預設情況下，設計器要求 Kozak ≥ 0.7 且 3'UTR ≥ 0.5。若沒有任何組合同時
滿足兩個閾值，則退回選取整體排名最高的組合。

```python
# 對 Kozak 嚴格、對 3'UTR 寬鬆
cfg = UTRDesignConfig(
    cds="MVSKGEELFTGV",
    prefer_kozak=0.85,
    prefer_utr3=0.3,
)
```

## 選擇 CDS 後端

預設情況下，CDS 以多目標後端（v0.25.0 RNop 知識注入式損失）進行優化。
`basic` 後端亦可使用：

```python
cfg = UTRDesignConfig(cds="MVSKGEELFTGV", backend="basic")
```

多目標後端之所以被推薦，是因為它在單一加權損失中同時考量
CAI + GC + CpG + 稀有連續 + 結構代理，給出比 CAI+GC 更細緻的 CDS 選擇。

## 輸出

結果是一個 `UTRDesignResult` 資料類別，內含組裝最終構件所需的全部資訊：

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

使用 `result.utr5 + result.cds_dna + result.utr3` 作為完整 mRNA 構件中
「CDS 兩側接 UTR」的部分。

## 未優化的部分

- **5'UTR / 3'UTR 的序列內容**：我們自手動策展的函式庫（4 + 3 = 7 個變體）
  中抽樣。若要做序列層級的 UTR 優化（這是另一個問題——在可微分
  表現量模型上做梯度式設計），請使用專用工具如 UTR-Function 或
  5UTR-Tailor。函式庫做法對不確定性更誠實，也能避免過擬合。
- **跨 UTR 配對的 CDS 逐位置選擇**：CDS 只優化一次，之後在 12 種 UTR 配對
  中重複使用。這在 CDS 與 UTR 選擇互相獨立的意義下是正確的，但也代表
  我們並未對 (CDS, UTR) 做聯合搜尋。實務上主導訊號來自 CDS 本身，
  所以這樣做沒有問題。
- **Poly-A 尾長度**：完整的 poly-A 組裝請見 `construct` 工具。
  `utr-design` 專注於 UTR × CDS 的耦合。