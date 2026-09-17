# 後端

七個工具，各自帶有具型別的 `Protocol` 配接器與純標準函式庫 mock 後備。
重型依賴依需求透過 `pip install` extras 啟用——CI 在沒有它們的情況下執行。

## 後端矩陣

| Backend | 工具 | 功能 | 安裝 |
|---|---|---|---|
| `mock` | 全部 | 純標準函式庫樁。永遠可用。 | — |
| `openai` | `neoantigen`、`trial` | 呼叫 OpenAI Chat Completions（或任何 OpenAI 相容端點） | `OPENAI_API_KEY=...` |
| `mhcflurry` | `neoantigen` | 全等位基因 MHC-I 結合親和力 | `pip install -e ".[neoantigen-mhcflurry]"` |
| `MedCPT` | `neoantigen`、`trial` | 生物醫學密集檢索（PubMed 對比學習） | `pip install -e ".[neoantigen-medcpt]"` 或 `[trial-medcpt]` |
| `scGPT` | `scrna` | 單細胞基礎模型嵌入 | `pip install -e ".[scrna]"` |
| `AlphaMissense` | `variant_scorer` | 透過 71M 變異 TSV 預測致變性 | 獨立 TSV（CC BY-NC-SA） |
| `AlphaGenome Atlas` | `variant-regulatory`、`variant_scorer`、`scrna` | 全部 9 億個可能 SNV 的調控變異影響（AVI） | `pip install -e \".[variant-alphagenome]\"` + `ALPHAGENOME_API_KEY` |
| `PhyloP46way` | `variant_scorer`、`scrna` | 每個 hg38 鹼基的演化保守度（UCSC 46-way 胎盤哺乳類比對） | （使用標準函式庫 urllib；無額外依賴） |
| `ESM2` | `neoantigen` | 用於免疫原性的凍結蛋白質 LM 嵌入 | `pip install -e ".[protein-lm]"` |
| `RiboDecode`（真實） | `codon` | 翻譯 × MFE 聯合密碼子最佳化 | `pip install ribodecode-1.3.0-py3-none-any.whl` |
| `RiboDecode`（啟發式） | `codon` | 純標準函式庫 RiboDecode 風格山丘爬山 | — |
| `LinearDesign` | `codon` | 翻譯 × MFE 聯合 DP | —（純標準函式庫） |
| `STModule`（真實） | `spatial` | 從 SRT 資料識別組織模組 | `pip install STModule` + Rscript 在 `$PATH` |
| `STModule`（mock） | `spatial` | 各平台基因宇宙的純標準函式庫樁 | — |
| `TrialGPT` | `trial` | 每條標準 LLM 資格配對 | （使用 `openai` 後端） |
| `Sim-ICL` | `trial` | 透過 TF-IDF 餘弦挑選 top-K 示範 | —（純標準函式庫） |

## Protocol 配接器如何運作

每個真實模型後端皆封裝於一個 `runtime_checkable` Protocol 之中：

```python
@runtime_checkable
class TranslationPredictor(Protocol):
    def predict(self, cds: str, env: str = "HEK293T",
                custom_env_csv: Path | None = None) -> TranslationPrediction: ...
```

本工具組同時提供真實配接器（例如 subprocess 呼叫 `pred-translation`
的 `TranslationModelCLIAdapter`）與 mock（使用 CAI 衍生分數的
`MockTranslationPredictor`）。消費者看見相同的 Protocol；後端選擇器
依 `$PATH` 與已安裝依賴來選擇真實或 mock。

## 後端選擇器

| 工具 | 選擇器 | 真實或 mock 分派 |
|---|---|---|
| `codon` | （CLI 旗標 `--backend`） | 若 `ribo-decode` 在 `$PATH` 則 `ribodecode-real`，否則 `ribodecode`（啟發式） |
| `neoantigen` | `select_translation_predictor` / `select_codon_optimizer` | 已安裝 mhcflurry；否則若設定 OpenAI key 則為 OpenAI；否則 mock |
| `trial` | （CLI 旗標 `--backend`） | 若已設定 OpenAI key 則為 OpenAI；否則 mock |
| `scrna` | `embed_with_foundation_model` | 已安裝 scGPT；否則身分識別（無操作後備） |
| `manufacture` | — | 全標準函式庫（無 LLM） |
| `lnp` | — | 全標準函式庫（無 LLM） |
| `spatial` | `select_spatial_module_backend` | `Rscript` 在 `$PATH`；否則 mock |

## 自動偵測

當 `--backend auto`（預設值）時，工具依下列順序挑選最強的可用後端：

1. **重型上游二進位**（例如 STModule 用 `Rscript`、RiboDecode 用
   `pred-translation`）——最佳準確度，需要使用者設定。
2. **已安裝的 Python 套件**（mhcflurry、transformers、OpenAI）
   ——良好的準確度，選擇性啟用。
3. **Mock**——永遠可用，確定性，~0 ms，僅標準函式庫。

這讓本機開發輕鬆，並讓正式部署透過環境變數或 CLI 旗標覆寫。
