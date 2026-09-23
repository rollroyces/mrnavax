# `construct` — 全 mRNA 構件組裝

從蛋白質胺基酸序列組裝完整的治療性 mRNA 構件（5'UTR + CDS + 3'UTR + 多聚腺苷酸化訊號 + poly-A 尾）。純標準庫——無重型依賴。

構件設計器是密碼子優化器的自然延伸：它不只是單獨優化 CDS，而是產生 mRNA 療法實際遞送給患者的完整構件。它組合了現有的基元：

* **5'UTR** — 帶強 Kozak context 的人類共識序列（`GGGCGACGCGGTGGCGGCCGCTCATGG`）。
* **CDS** — 從蛋白質 AA 序列反向翻譯，再透過 [codon](../tools/codon.md) 後端之一進行優化（`basic` / `ribodecode` / `lineardesign` / `ribodecode-real` / **`multi-objective` ← 自 v0.25.0 起為預設**）。
* **3'UTR** — 帶兩個 ARE 穩定元素的人類共識序列。
* **多聚腺苷酸化訊號** — `AAUAAA`（DNA：`AATAAA`）。
* **Poly-A 尾** — 預設 120 nt（可透過 `--poly-a-length` 設定）。

## 用法

```bash
# 從蛋白質 AA 序列組裝完整 mRNA 構件
mrnavax construct --sequence "MVSKGEELFTGV"

# 使用特定的 CDS 優化後端
mrnavax construct --sequence my_protein.fasta \
    --backend lineardesign

# 調整 poly-A 尾長度（治療性 mRNA 典型值：100-150 nt）
mrnavax construct --sequence my_protein.fasta \
    --poly-a-length 100

# 將完整報告儲存為 JSON 檔案
mrnavax construct --sequence my_protein.fasta \
    --out my_construct.json
```

## 為何重要

已發表的治療性 mRNA 設計（Moderna mRNA-1273、BioNTech BNT162b2）共享同一個共同架構：5'UTR（Kozak context）+ CDS + 3'UTR（ARE 元素）+ 多聚腺苷酸化訊號 + poly-A 尾。在 v0.26.0 之前，mrnavax 將各個元件單獨暴露出來（`codon` 處理 CDS、`manufacture` 處理 QC 檢查）。自 v0.26.0 起，`construct` 工具能以一次 CLI 呼叫組裝完整的治療性 mRNA 架構，並預設對 CDS 套用 SOTA 多目標優化器。

## 參考文獻

5'/3' UTR 模板遵循已發表治療性 mRNA 設計（Moderna mRNA-1273、BioNTech BNT162b2）的人類共識序列。CDS 優化後端與 `mrnavax codon` 的五個後端相同。