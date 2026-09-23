# `construct` — 全 mRNA 构建组装

从蛋白质氨基酸序列组装完整的治疗性 mRNA 构建（5'UTR + CDS + 3'UTR + 多聚腺苷酸化信号 + poly-A 尾）。纯标准库——无重型依赖。

构建设计器是密码子优化器的自然延伸：它不只是单独优化 CDS，而是产生 mRNA 疗法实际递送给患者的完整构建。它组合了现有的基元：

* **5'UTR** — 带强 Kozak context 的人类共识序列（`GGGCGACGCGGTGGCGGCCGCTCATGG`）。
* **CDS** — 从蛋白质 AA 序列反向翻译，再通过 [codon](../tools/codon.md) 后端之一进行优化（`basic` / `ribodecode` / `lineardesign` / `ribodecode-real` / **`multi-objective` ← 自 v0.25.0 起为默认**）。
* **3'UTR** — 带两个 ARE 稳定元素的人类共识序列。
* **多聚腺苷酸化信号** — `AAUAAA`（DNA：`AATAAA`）。
* **Poly-A 尾** — 默认 120 nt（可通过 `--poly-a-length` 设置）。

## 用法

```bash
# 从蛋白质 AA 序列组装完整 mRNA 构建
mrnavax construct --sequence "MVSKGEELFTGV"

# 使用特定的 CDS 优化后端
mrnavax construct --sequence my_protein.fasta \
    --backend lineardesign

# 调整 poly-A 尾长度（治疗性 mRNA 典型值：100-150 nt）
mrnavax construct --sequence my_protein.fasta \
    --poly-a-length 100

# 将完整报告保存为 JSON 文件
mrnavax construct --sequence my_protein.fasta \
    --out my_construct.json
```

## 为何重要

已发表的治疗性 mRNA 设计（Moderna mRNA-1273、BioNTech BNT162b2）共享同一个共同架构：5'UTR（Kozak context）+ CDS + 3'UTR（ARE 元素）+ 多聚腺苷酸化信号 + poly-A 尾。在 v0.26.0 之前，mrnavax 将各个元件单独暴露出来（`codon` 处理 CDS、`manufacture` 处理 QC 检查）。自 v0.26.0 起，`construct` 工具能以一次 CLI 调用组装完整的治疗性 mRNA 架构，并默认对 CDS 应用 SOTA 多目标优化器。

## 参考文献

5'/3' UTR 模板遵循已发表治疗性 mRNA 设计（Moderna mRNA-1273、BioNTech BNT162b2）的人类共识序列。CDS 优化后端与 `mrnavax codon` 的五个后端相同。