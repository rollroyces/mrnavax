# 后端

七个工具，各自带有具类型的 `Protocol` 适配器与纯标准库 mock 后备。
重型依赖依需求通过 `pip install` extras 启用——CI 在没有它们的情况下运行。

## 后端矩阵

| Backend | 工具 | 功能 | 安装 |
|---|---|---|---|
| `mock` | 全部 | 纯标准库桩。永远可用。 | — |
| `openai` | `neoantigen`、`trial` | 调用 OpenAI Chat Completions（或任何 OpenAI 兼容端点） | `OPENAI_API_KEY=...` |
| `mhcflurry` | `neoantigen` | 全等位基因 MHC-I 结合亲和力 | `pip install -e ".[neoantigen-mhcflurry]"` |
| `MedCPT` | `neoantigen`、`trial` | 生物医学密集检索（PubMed 对比学习） | `pip install -e ".[neoantigen-medcpt]"` 或 `[trial-medcpt]` |
| `scGPT` | `scrna` | 单细胞基础模型嵌入 | `pip install -e ".[scrna]"` |
| `AlphaMissense` | `variant_scorer` | 通过 71M 变异 TSV 预测致变性 | 独立 TSV（CC BY-NC-SA） |
| `AlphaGenome Atlas` | `variant-regulatory`、`variant_scorer`、`scrna` | 全部 90 亿个可能 SNV 的调控变异影响（AVI） | `pip install -e \".[variant-alphagenome]\"` + `ALPHAGENOME_API_KEY` |
| `ESM2` | `neoantigen` | 用于免疫原性的冻结蛋白质 LM 嵌入 | `pip install -e ".[protein-lm]"` |
| `RiboDecode`（真实） | `codon` | 翻译 × MFE 联合密码子优化 | `pip install ribodecode-1.3.0-py3-none-any.whl` |
| `RiboDecode`（启发式） | `codon` | 纯标准库 RiboDecode 风格爬山 | — |
| `LinearDesign` | `codon` | 翻译 × MFE 联合 DP | —（纯标准库） |
| `STModule`（真实） | `spatial` | 从 SRT 数据识别组织模块 | `pip install STModule` + Rscript 在 `$PATH` |
| `STModule`（mock） | `spatial` | 各平台基因宇宙的纯标准库桩 | — |
| `TrialGPT` | `trial` | 每条标准 LLM 资格配对 | （使用 `openai` 后端） |
| `Sim-ICL` | `trial` | 通过 TF-IDF 余弦挑选 top-K 示范 | —（纯标准库） |

## Protocol 适配器如何工作

每个真实模型后端皆封装于一个 `runtime_checkable` Protocol 之中：

```python
@runtime_checkable
class TranslationPredictor(Protocol):
    def predict(self, cds: str, env: str = "HEK293T",
                custom_env_csv: Path | None = None) -> TranslationPrediction: ...
```

本工具包同时提供真实适配器（例如 subprocess 调用 `pred-translation`
的 `TranslationModelCLIAdapter`）与 mock（使用 CAI 衍生分数的
`MockTranslationPredictor`）。消费者看见相同的 Protocol；后端选择器
依 `$PATH` 与已安装依赖来选择真实或 mock。

## 后端选择器

| 工具 | 选择器 | 真实或 mock 分派 |
|---|---|---|
| `codon` | （CLI 标志 `--backend`） | 若 `ribo-decode` 在 `$PATH` 则 `ribodecode-real`，否则 `ribodecode`（启发式） |
| `neoantigen` | `select_translation_predictor` / `select_codon_optimizer` | 已安装 mhcflurry；否则若设置 OpenAI key 则为 OpenAI；否则 mock |
| `trial` | （CLI 标志 `--backend`） | 若已设置 OpenAI key 则为 OpenAI；否则 mock |
| `scrna` | `embed_with_foundation_model` | 已安装 scGPT；否则身份识别（无操作后备） |
| `manufacture` | — | 全标准库（无 LLM） |
| `lnp` | — | 全标准库（无 LLM） |
| `spatial` | `select_spatial_module_backend` | `Rscript` 在 `$PATH`；否则 mock |

## 自动检测

当 `--backend auto`（默认值）时，工具依下列顺序挑选最强的可用后端：

1. **重型上游二进制**（例如 STModule 用 `Rscript`、RiboDecode 用
   `pred-translation`）——最佳准确度，需要用户设置。
2. **已安装的 Python 包**（mhcflurry、transformers、OpenAI）
   ——良好的准确度，选择性启用。
3. **Mock**——永远可用，确定性，~0 ms，仅标准库。

这让本地开发轻松，并让正式部署通过环境变量或 CLI 标志覆盖。
