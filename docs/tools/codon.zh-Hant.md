# `codon` — 序列分析與優化

計算所有現代 mRNA 設計模型（CodonBERT、RiboDecode、LinearDesign、mRNABERT）所共同仰賴的標準密碼子使用特徵。

## 用法

```bash
mrnavax codon --sequence mrnavax/examples/cas9.fasta
mrnavax codon --sequence mrnavax/examples/cas9.fasta --optimize
```

## 輸出 schema

```json
{
  "n_codons": 210,
  "cai": 0.6954,                        // Codon Adaptation Index (0–1)
  "gc_percent": 37.3,                   // 整體 GC%
  "rare_codon_fraction": 0.0429,        // 使用頻率低於 0.10 的密碼子比例
  "cpg_obs_exp": 1.1331,                // CpG O/E 比值——先天免疫之代理指標
  "most_common_codons": [["GAT", 12], ...],
  "rare_codons": ["CTA", "TTA"],
  "gc_window_stddev": 4.75               // 滑動 GC 標準差（轉譯速度之代理指標）
}
```

## `--optimize`

執行一輪貪婪式的同義密碼子替換，在維持 GC% 於 45–60% 區間的前提下，極大化各密碼子的使用頻率。

對一條在人類細胞表現的細菌基因，預期 CAI 約可提升 0.2（例如 Cas9 範例由 0.7 → 0.93）。

這是任何現代密碼子模型都應該超越的**經典基準**。刻意保持簡單，以便您能將更精密的模型（CodonBERT、RiboDecode）以相同的輸入輸出形狀接入並進行比較。

## 多目標優化（知識注入式損失）

`multi-objective` 後端實現了來自 arXiv:2505.23862（RNop，*mRNA Design and Optimization with Deep Knowledge-Infused Approach*，2026 年 8 月）的 SOTA 模式：將每個生物學目標分解為各自的評分，使使用者能看到**每個**同義密碼子替換是由**哪一個**分量驅動的。這正是頂級專案（RNop / LinearDesign / RiboDecode）共同匯聚的「不可能三角形」（保真度、多目標、效率）解法：保真度（僅同義替換）、多目標明確（加權損失向量）、效率以純 Python 標準庫達成。

五個歸一化分量：

| 分量 | 範圍 | 衡量內容 |
|---|---|---|
| `cai` | (0, 1] | Codon Adaptation Index — 轉譯效率 |
| `gc_score` | [0, 1] | GC% 是否落入目標區間 — 轉錄穩定性之代理指標 |
| `cpg_score` | [0, 1] | CpG Obs/Exp 規避 — 先天免疫之代理指標 |
| `rare_run_penalty` | [0, 1] | 最長稀有密碼子連續長度 — 核糖體停滯之代理指標 |
| `structure_proxy` | [0, 1] | GC 滑動窗口標準差的一致性 — 局部結構之代理指標 |

整體分數為各分量的加權總和（`rare_run_penalty` 以減項處理），夾至 `[0, 1]`。預設權重遵循 SOTA 典型模式：CAI 為主導訊號（~40%），其次為結構代理（~25%），再來是 GC 與 CpG，稀有密碼子連續長度計為軟性懲罰。

```bash
mrnavax codon --sequence gfp.fasta \
    --optimize --backend multi-objective
```

輸出包含 `before` 與 `after` 各分量的分解，加上 `overall_before`、`overall_after`、`improvement` 與 `n_changes`，使各生物軸的貢獻完全可審計。

**自訂損失向量**：`MultiObjectiveConfig` 暴露全部 5 個權重（以及 GC 區間、稀有閾值、結構窗口）。在 Python 中可自行配置：

```python
from mrnavax.codon_multi_objective import (
    MultiObjectiveConfig, multi_objective_optimize,
)

cfg = MultiObjectiveConfig(
    cai_weight=0.30,        # 略降 CAI 比重，提高結構權重
    gc_weight=0.20,
    cpg_weight=0.15,
    rare_run_weight=0.25,   # 強調避免稀有密碼子連續
    structure_weight=0.10,
    target_gc_min=50.0,     # 更緊的 GC 區間
    target_gc_max=58.0,
)
result = multi_objective_optimize(cds, cfg)
print(result.improvement, result.n_changes)
```

## RNop 參考文獻

Gong, Z., Jiang, Z., Gao, W., Wang, Y., Cai, Z., Deng, Z., Ma, L.
(2026). *mRNA Design and Optimization with Deep Knowledge-Infused
Approach.* arXiv:2505.23862v2。「不可能三角形」的論述（保真度、多目標、效率）與知識注入式損失分解為 `multi-objective` 後端分量分解的概念基礎。
