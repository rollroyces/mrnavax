# 贡献指南

欢迎提交 Pull Request。

## 限制

- **纯标准库核心。** 在默认安装（`pip install -e .`）中随附的
  工具必须仅依赖 Python 标准库。新依赖须置于 `pyproject.toml`
  的 `optional-dependencies` extra 之后。
- **确定性。** 工具对相同输入必须产生相同输出（默认输出不含时间戳）。
  Mock 后端永远是确定性的。
- **结构验证。** 由 LLM 支持的工具必须验证响应形状，若 LLM 返回
  格式错误的 JSON 则退回启发式方法。
- **精简 CLI 界面。** 每个工具的 `_run_cli` 应接受 `--out`，未设置时
  输出至 stdout。
- **基于 Protocol 的可插拔性。** 重型模型集成（RiboDecode、STModule、
  ESM2、scGPT、mhcflurry）必须置于 `runtime_checkable` Protocol 或抽象
  适配器类之后，通过可选 extras 安装（例如
  `pip install mrnavax[sota]`）。参考
  `mrnavax/codon_ribodecode_adapter.py`、
  `mrnavax/spatial_module_adapter.py`、
  `mrnavax/protein_lm_adapter.py`。
- **新功能采用 TDD。** 每个新的公开函数必须在实现前先写测试
  （见下方 *开发* 段落）。

## 开发

```bash
git clone https://github.com/rollroyces/mrnavax.git
cd mrnavax
pip install -e ".[dev,llm]"

# 运行 25 项后端完整性检查
python -m mrnavax.backends --check-all

# 运行单元测试套件（167 个测试）
python -m unittest discover tests

# 烟雾测试
bash scripts/smoke.sh

# 本地文档
pip install -e ".[docs]"
mkdocs serve
```

## Pull Request 流程

1. 非显而易见的修改请先开 issue。
2. 新功能遵循 **TDD 循环**（依 `test-driven-development` 技能）：
   - **RED**：编写一个失败的测试，练习期望的 API。运行它并确认
     因正确的原因失败。
   - **GREEN**：编写最小实现使其通过。
   - **REFACTOR**：清理重复、命名、辅助函数。
3. 包装第三方模型时遵循 **api-integration-verify** 纪律：先探查
   上游文档／GitHub，从真实 URL 验证参数名称与返回形状，再设计
   Protocol 适配器。
4. 在 `backends.py` 中为任何新后端新增 `register()` 条目。
5. CI 必须在 Python 3.11–3.14 上保持绿灯才能合并。
6. 当 `bash scripts/smoke.sh` +
   `python -m mrnavax.backends --check-all` +
   `python -m unittest discover tests` 全部通过后，使用 `[verified]`
   commit 信息进行 squash-merge。

## 新增工具

每个新工具应随附：

1. 输入与输出的具类型 `dataclass`（frozen，于构造时验证）。
2. 后端接口的 `runtime_checkable` Protocol。
3. 通过 subprocess / 延迟加载接入上游模型的真实适配器。
4. 满足相同 Protocol 的纯标准库 mock。
5. `backends.py` 中的 `register()` 条目供 CI 完整性使用。
6. `tests/` 中遵循严格 TDD 的测试。
7. 一份 `docs/tools/<name>.{md,zh-Hant.md,zh-Hans.md}` 页面，描述
   工具、CLI 用法、后端矩阵与参考论文。

## 重新整理记录的 Atlas fixture

套件中的 fixture 位于 `tests/fixtures/alphagenome_atlas_sample.json`，
会锁定解析器形状以防上游变更。若要使用真实撷取响应重新整理（取得
`ALPHAGENOME_API_KEY` 后）：

```bash
# 1. 导出您的 API 密钥（https://deepmind.google.com/science/alphagenome）
export ALPHAGENOME_API_KEY=***

# 2. 手动触发 GitHub Actions 工作流程：
#    Actions → Atlas live integration → Run workflow
#    Inputs: mode=record, leave baseline/tolerance as default
#    这会执行真实的 Atlas 调用并将响应提交为
#    tests/fixtures/alphagenome_atlas_live_<timestamp>.json。
#
# 3. 审阅 PR，将撷取响应复制到套件 fixture（或一并提交），并更新
#    任何整合测试中的预期评分断言。
#
# 4. 本机也可以执行：
python -m mrnavax.live_atlas_integration --mode record \
    --output tests/fixtures/alphagenome_atlas_live.json
```

每周排程（`.github/workflows/atlas_integration.yml`）以 `--mode regression`
执行同一个测试工具，并以 ±5% 容差将实时评分与基准比较——当评分
漂移超出容差时工作流程会失败。这能即时抓到上游 Atlas API 的损坏
（记录的 fixture 则会有 7 天的延迟）。

## 考虑过但未构建

以下功能已经过 **试探且刻意延后**。每个条目记录了问题、裁决，
以及日后应该重新审视的条件。未来的贡献者在提出类似功能前应该
先阅读本节——每个「明显的」版本都已被考虑过并附理由拒绝。

### `docs/live-atlas.md` — 自动更新的每周结果页（v0.24.0 试探）

**问题：** 每周 Atlas 排程是否也应该将结果发布到公开文档站，
作为「最近一次实时执行 + 历史记录」页面？

**裁决：** PARTIAL → **延后**（最小可行变更约 150 行代码——
工作流程 + 脚本 + 页面；在 API 密钥设置完成前，成本高于效益）。

**反面案例（延后原因）：**

1. 此存储库尚未设置 API 密钥，因此短期内页面会向全世界显示
   自己的停用状态。我们并未提供可见的价值。
2. 7 天的执行周期加上 v0.18.0 记录的 fixture 已经能抓到 Atlas
   的损坏。这个页面只在「有东西坏掉时」增加可见性——但在大多
   数正常的周次，表格只会是 ~52 列几乎相同的「✅ within ±0.05」。
3. 成本：新增一支脚本、一个工作流程工作、持续维护，以及在排程
   上增加 `contents: write` 权限。

**重新审视时机：** （a）已设置 `ALPHAGENOME_API_KEY` GitHub
   secret，且（b）我们希望有一个用户可见的「活跃 vs 弃用」
   信号。完整的试探设计已保存于
   `.considerations/live-atlas-docs-page.md`（位于 repo 根目录；
   或搜索 git 记录中的 v0.24.0 试探 commit）。

### 更多条目将陆续新增

凡是通过此模式拒绝的功能，请在关闭该试探的同一个 commit 中，
将裁决记录到这里。未来的读者应该可以 grep CONTRIBUTING 中的
「考虑过但未构建」关键字，看到所有曾被试探但未出货的功能请求
完整历史。

## 许可

贡献代码即代表您同意自己的贡献将依本项目的双重许可条款许可
（AGPL-3.0-or-later 用于开源使用；商业许可可来信索取）。
