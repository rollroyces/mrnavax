# `predict` — 即時 AlphaGenome Atlas 調控變異預測

mrnavax 工具組中的**第 11 個工具**。呼叫真正的 AlphaGenome
Atlas API，對一個或多個變異的調控影響（AVI 分數）進行評分。
這是工具組中**預設會發起即時網路呼叫**的唯一工具；其餘全部在本地執行。

## 為何要獨立成一個 CLI 工具？

v0.25.1 版本加入了對 `ALPHAGENOME_API_KEY` 的輔助檔案自動偵測。
加上這個之後，透過即時 Atlas API 對一個變異進行端到端評分便成了一行可解。
`predict` 工具將這一行包裝起來，提供：

- 單變異與 CSV 批次模式
- 自動 API 金鑰解析（環境變數 / 輔助檔案）
- 健全的子行程錯誤處理
- 適合下游管線的 JSON 輸出

## 快速上手

### 單一變異

```bash
mrnavax predict --chrom chr7 --pos 140753336 --ref T --alt A --name BRAF_V600E
```

### 從 CSV 批次

```csv
name,chrom,pos,ref,alt
BRAF_V600E,chr7,140753336,T,A
KRAS_G12D,chr12,25245350,C,T
```

```bash
mrnavax predict --input variants.csv
```

### 輸出

```json
{
  "n_variants": 1,
  "results": [
    {
      "name": "BRAF_V600E",
      "chrom": "chr7",
      "pos": 140753336,
      "ref": "T",
      "alt": "A",
      "score": 1.0,
      "classification": "high",
      "is_coding": true,
      "n_scorers": 19
    }
  ]
}
```

## 輸出欄位

|| 欄位 | 型別 | 意義 |
||---|---|---|
|| `score` | float in [0, 1] | 19 個 Atlas 評分器的最大絕對差值，裁剪至 [0, 1]。值越高代表調控影響越大。 |
|| `classification` | str in {"low", "moderate", "high"} | 標準 AVI 分箱（low < 0.34，moderate < 0.564，否則為 high）。 |
|| `is_coding` | bool | 若變異位置與 16kb 視窗內的 protein_coding 基因重疊則為 True。屬最佳近似；請視為提示而非精確註解。 |
|| `n_scorers` | int | 有回應的 Atlas 評分器數量（通常為 19）。 |

## API 金鑰解析

shim 按下列順序自動解析 API 金鑰：

1. `api_key=...` 明確指定（於 JSON payload 中——用於直接程式化呼叫）
2. `ALPHAGENOME_API_KEY` 環境變數（CI 路徑）
3. `~/projects/alphagenome-work/.alphagenome_key` 輔助檔案（本地開發）
4. `~/.alphagenome_key` 替代位置

若以上皆未解析成功，shim 會以 exit code 1 結束，並在 stderr
回傳 `"API key not configured"`。

## CLI vs. 程式化

程式化使用建議採 Python API：

```python
import json, subprocess
from pathlib import Path

shim = Path("/path/to/mrnavax/_shims/alphagenome_cli.py")
proc = subprocess.run(
    [sys.executable, str(shim)],
    input=json.dumps({"chrom": "chr7", "pos": 140753336, "ref": "T", "alt": "A"}),
    capture_output=True, text=True, timeout=180,
)
result = json.loads(proc.stdout)
print(result["score"], result["classification"])
```

`predict` CLI 不過是這個子行程呼叫的薄包裝。

## 效能

- 單一變異：約 3–5 秒（API 暖機後）。冷啟動首次呼叫可能約 30 秒。
- 批次：循序執行；每個變異約 3–5 秒。內建未做速率限制處理；
  若需處理上千個變異，請直接透過 Python API 並自行控制並行。

## 成本 / 配額

Atlas API 對學術用途免費（預覽階段）；執行大批次前請查閱
最新的 Google AI Studio Atlas 條款。本工具組對配額不做任何聲明；
shim 會原樣傳遞上游 API 回傳的結果。

## `predict` *不是*什麼

- **不是致病性預測器**：AVI 分數衡量的是調控影響，並非疾病成因。
  若要做診斷性宣稱，請搭配 AlphaMissense（編碼變異）+ 臨床遺傳學詮釋。
- **不是結構影響預測器**：蛋白質結構 / 功能效應請使用 AlphaMissense。
- **不是多基因風險分數**：Atlas 報告的是單一變異的影響；
  PRS 需要彙整多個變異。

## 參考文獻

- Avsec et al., "AlphaGenome: AI for genomics at scale", *Nature*
  2026（DOI 待出版後補上）。
- https://deepmind.google.com/science/alphagenome — 官方 API 文件。