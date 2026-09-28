# `predict` — 实时 AlphaGenome Atlas 调控变异预测

mrnavax 工具组中的**第 11 个工具**。调用真正的 AlphaGenome
Atlas API，对一个或多个变异的调控影响（AVI 分数）进行评分。
这是工具组中**默认会发起实时网络调用**的唯一工具；其余全部在本地执行。

## 为何要独立成一个 CLI 工具？

v0.25.1 版本加入了 `ALPHAGENOME_API_KEY` 的辅助文件自动侦测。
加上这个之后，通过实时 Atlas API 对一个变异进行端到端评分便成了一行可解。
`predict` 工具将这一行包装起来，提供：

- 单变异与 CSV 批次模式
- 自动 API 密钥解析（环境变量 / 辅助文件）
- 健壮的子进程错误处理
- 适合下游管道的 JSON 输出

## 快速上手

### 单个变异

```bash
mrnavax predict --chrom chr7 --pos 140753336 --ref T --alt A --name BRAF_V600E
```

### 从 CSV 批次

```csv
name,chrom,pos,ref,alt
BRAF_V600E,chr7,140753336,T,A
KRAS_G12D,chr12,25245350,C,T
```

```bash
mrnavax predict --input variants.csv
```

### 输出

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

## 输出字段

|| 字段 | 类型 | 含义 |
||---|---|---|
|| `score` | float in [0, 1] | 19 个 Atlas 评分器的最大绝对差值，裁剪至 [0, 1]。值越高代表调控影响越大。 |
|| `classification` | str in {"low", "moderate", "high"} | 标准 AVI 分箱（low < 0.34，moderate < 0.564，否则为 high）。 |
|| `is_coding` | bool | 若变异位置与 16kb 窗口内的 protein_coding 基因重叠则为 True。属最佳近似；请视为提示而非精确注释。 |
|| `n_scorers` | int | 有响应的 Atlas 评分器数量（通常为 19）。 |

## API 密钥解析

shim 按以下顺序自动解析 API 密钥：

1. `api_key=...` 显式指定（于 JSON payload 中——用于直接编程化调用）
2. `ALPHAGENOME_API_KEY` 环境变量（CI 路径）
3. `~/projects/alphagenome-work/.alphagenome_key` 辅助文件（本地开发）
4. `~/.alphagenome_key` 备用位置

若以上皆未解析成功，shim 会以 exit code 1 退出，并在 stderr
返回 `"API key not configured"`。

## CLI vs. 编程化

编程化使用建议采用 Python API：

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

`predict` CLI 不过是这个子进程调用的薄包装。

## 性能

- 单个变异：约 3–5 秒（API 暖机后）。冷启动首次调用可能约 30 秒。
- 批次：顺序执行；每个变异约 3–5 秒。内置未做速率限制处理；
  若需处理上千个变异，请直接通过 Python API 并自行控制并发。

## 成本 / 配额

Atlas API 对学术用途免费（预览阶段）；执行大批次前请查阅
最新的 Google AI Studio Atlas 条款。本工具组对配额不作任何声明；
shim 会原样传递上游 API 返回的结果。

## `predict` *不是*什么

- **不是致病性预测器**：AVI 分数衡量的是调控影响，并非疾病成因。
  若要做诊断性宣称，请搭配 AlphaMissense（编码变异）+ 临床遗传学诠释。
- **不是结构影响预测器**：蛋白质结构 / 功能效应请使用 AlphaMissense。
- **不是多基因风险分数**：Atlas 报告的是单个变异的影响；
  PRS 需要聚合多个变异。

## 参考文献

- Avsec et al., "AlphaGenome: AI for genomics at scale", *Nature*
  2026（DOI 待出版后补上）。
- https://deepmind.google.com/science/alphagenome — 官方 API 文档。