# Atlas API 密钥启用

工具包默认带有真实数据的 AlphaGenome Atlas 集成，但在实时每周定时任务可以拉取真实数据之前，需要先配置 Google AI Studio API 密钥。本页涵盖**三种启用路径**：

## 为什么要三种路径？

密钥最多可以存在三个位置：

  1. `ALPHAGENOME_API_KEY` **环境变量** — CI 首选
  2. **`~/projects/alphagenome-work/.alphagenome_key`** — 本地开发首选（自动检测；v0.25.1 新增）
  3. **`~/.alphagenome_key`** — 备用 home 位置（v0.25.1 新增）

加载器按此顺序解析：显式的 `api_key=...` 参数 > `ALPHAGENOME_API_KEY` 环境变量 > helper file `~/projects/alphagenome-work/.alphagenome_key` > `~/.alphagenome_key`。

## 路径 1：本地开发（推荐）

如果你 clone 了此 repo 并希望 Atlas 调用从笔记本上「直接可用」，请将 39 字符的 Google AI Studio 密钥（以 `AIza` 开头）写入：

```bash
mkdir -p ~/projects/alphagenome-work
echo -n 'YOUR_39_CHAR_KEY_HERE' > ~/projects/alphagenome-work/.alphagenome_key
chmod 600 ~/projects/alphagenome-work/.alphagenome_key
```

下次运行 `score_variant()` 或 `live_atlas_integration.py` 时，加载器会自动检测此文件。

## 路径 2：GitHub Actions 每周定时任务（CI 推荐）

要启用每周的 Atlas 集成测试：

  1. 从 https://aistudio.google.com/apikey 获取 Google AI Studio API 密钥（提供免费套餐）。
  2. 打开 https://github.com/rollroyces/mrnavax/settings/secrets/actions/new
  3. 名称：`ALPHAGENOME_API_KEY`
  4. 值：粘贴 39 字符的密钥
  5. 点击「Add secret」

下次定时 cron（每周一 06:00 UTC）或手动 `gh workflow run atlas_integration.yml` 将会调用实时 Atlas API，并将结果与 `tests/fixtures/alphagenome_atlas_live.json` 中的 fixture 进行比对。如果上游 API 契约变更，集成测试将在一周内捕获。

## 路径 3：CI / 脚本（环境变量）

对于临时环境（Docker、没有 helper file 的 CI runner）：

```bash
export ALPHAGENOME_API_KEY='YOUR_39_CHAR_KEY_HERE'
```

## 验证启用

配置任何路径后，运行验证脚本以确认端到端 Atlas 连接：

```bash
python scripts/check_atlas_activation.py
```

预期输出（hg38 的 BRAF V600E）：

```
[helper-file] OK: /Users/hermes/projects/alphagenome-work/.alphagenome_key (40 bytes, prefix matches Google AI key: True)
[package] alphagenome installed
[key] OK: resolved 39-char key (prefix matches Google AI: True)
[live-call] score=1.000, classification='high', is_coding=True, n_scorers=19
[fixture] recorded=1.000/'high', live=1.000/'high', matches=True
[result] PASS: live Atlas matches recorded fixture
[result] Activation verified — weekly Atlas cron will work
```

如果实时调用失败，脚本将打印上游错误消息并以非零状态退出。

## 「已启用」的定义

Atlas API 密钥「已启用」意味着**全部三项**都为真：

  1. ✅ 密钥能通过上述三种路径之一成功解析
  2. ✅ alphagenome 软件包已安装（`pip install 'mrnavax[variant-alphagenome]'`）
  3. ✅ 实时 BRAF V600E 调用返回 `score≈1.0`、`classification='high'`、`is_coding=True`（与记录的 fixture 一致）

当三项都为真时，GitHub Actions 每周定时任务将执行实时 Atlas 调用，而非以「ALPHAGENOME_API_KEY not configured」通知跳过。

## 隐私

API 密钥永远不会被记录、回显或由工具包写入磁盘。验证脚本仅打印其长度与前缀匹配状态，而非密钥本身。