# 貢獻指南

歡迎提交 Pull Request。

## 限制

- **純標準函式庫核心。** 在預設安裝（`pip install -e .`）中隨附的
  工具必須僅依賴 Python 標準函式庫。新依賴須置於 `pyproject.toml`
  的 `optional-dependencies` extra 之後。
- **確定性。** 工具對相同輸入必須產生相同輸出（預設輸出不含時間戳）。
  Mock 後端永遠是確定性的。
- **結構驗證。** 由 LLM 支援的工具必須驗證回應形狀，若 LLM 傳回
  格式錯誤的 JSON 則退回啟發式方法。
- **精簡 CLI 介面。** 每個工具的 `_run_cli` 應接受 `--out`，未設定時
  輸出至 stdout。
- **基於 Protocol 的可插拔性。** 重型模型整合（RiboDecode、STModule、
  ESM2、scGPT、mhcflurry）必須置於 `runtime_checkable` Protocol 或抽象
  配接器類別之後，透過選用 extras 安裝（例如
  `pip install mrnavax[sota]`）。參考
  `mrnavax/codon_ribodecode_adapter.py`、
  `mrnavax/spatial_module_adapter.py`、
  `mrnavax/protein_lm_adapter.py`。
- **新功能採用 TDD。** 每個新的公開函式必須在實作前先寫測試
  （見下方 *開發* 段落）。

## 開發

```bash
git clone https://github.com/rollroyces/mrnavax.git
cd mrnavax
pip install -e ".[dev,llm]"

# 執行 25 項後端完整性檢查
python -m mrnavax.backends --check-all

# 執行單元測試套件（167 個測試）
python -m unittest discover tests

# 煙霧測試
bash scripts/smoke.sh

# 本機文件
pip install -e ".[docs]"
mkdocs serve
```

## Pull Request 流程

1. 非顯而易見的修改請先開 issue。
2. 新功能遵循 **TDD 循環**（依 `test-driven-development` 技能）：
   - **RED**：撰寫一個失敗的測試，練習期望的 API。執行它並確認
     因正確的原因失敗。
   - **GREEN**：撰寫最小實作使其通過。
   - **REFACTOR**：清理重複、命名、輔助函式。
3. 包裝第三方模型時遵循 **api-integration-verify** 紀律：先探查
   上游文件／GitHub，從真實 URL 驗證參數名稱與回傳形狀，再設計
   Protocol 配接器。
4. 在 `backends.py` 中為任何新後端新增 `register()` 條目。
5. CI 必須在 Python 3.11–3.14 上保持綠燈才能合併。
6. 當 `bash scripts/smoke.sh` +
   `python -m mrnavax.backends --check-all` +
   `python -m unittest discover tests` 全部通過後，使用 `[verified]`
   commit 訊息進行 squash-merge。

## 新增工具

每個新工具應隨附：

1. 輸入與輸出的具型別 `dataclass`（frozen，於建構時驗證）。
2. 後端介面的 `runtime_checkable` Protocol。
3. 透過 subprocess / 延遲載入接入上游模型的真實配接器。
4. 滿足相同 Protocol 的純標準函式庫 mock。
5. `backends.py` 中的 `register()` 條目供 CI 完整性使用。
6. `tests/` 中遵循嚴格 TDD 的測試。
7. 一份 `docs/tools/<name>.{md,zh-Hant.md,zh-Hans.md}` 頁面，描述
   工具、CLI 用法、後端矩陣與參考論文。

## 重新整理記錄的 Atlas fixture

套件中的 fixture 位於 `tests/fixtures/alphagenome_atlas_sample.json`，
會鎖定解析器形狀以防上游變更。若要使用真實擷取回應重新整理（取得
`ALPHAGENOME_API_KEY` 後）：

```bash
# 1. 匯出您的 API 金鑰（https://deepmind.google.com/science/alphagenome）
export ALPHAGENOME_API_KEY=***

# 2. 手動觸發 GitHub Actions 工作流程：
#    Actions → Atlas live integration → Run workflow
#    Inputs: mode=record, leave baseline/tolerance as default
#    這會執行真實的 Atlas 呼叫並將回應提交為
#    tests/fixtures/alphagenome_atlas_live_<timestamp>.json。
#
# 3. 審閱 PR，將擷取回應複製到套件 fixture（或一併提交），並更新
#    任何整合測試中的預期評分斷言。
#
# 4. 本機也可以執行：
python -m mrnavax.live_atlas_integration --mode record \
    --output tests/fixtures/alphagenome_atlas_live.json
```

每週排程（`.github/workflows/atlas_integration.yml`）以 `--mode regression`
執行同一個測試工具，並以 ±5% 容差將即時評分與基準比較——當評分
漂移超出容差時工作流程會失敗。這能即時抓到上游 Atlas API 的損壞
（記錄的 fixture 則會有 7 天的延遲）。

## 考慮過但未建構

以下功能已經過 **試探且刻意延後**。每個條目記錄了問題、裁決，
以及日後應該重新審視的條件。未來的貢獻者在提出類似功能前應該
先閱讀本節——每個「明顯的」版本都已被考慮過並附理由拒絕。

### `docs/live-atlas.md` — 自動更新的每週結果頁（v0.24.0 試探）

**問題：** 每週 Atlas 排程是否也應該將結果發布到公開文件站，
作為「最近一次即時執行 + 歷史紀錄」頁面？

**裁決：** PARTIAL → **延後**（最小可行變更約 150 行程式碼——
工作流程 + 腳本 + 頁面；在 API 金鑰設定完成前，成本高於效益）。

**反面案例（延後原因）：**

1. 此儲存庫尚未設定 API 金鑰，因此短期內頁面會向全世界顯示
   自己的停用狀態。我們並未提供可見的價值。
2. 7 天的執行週期加上 v0.18.0 記錄的 fixture 已經能抓到 Atlas
   的損壞。這個頁面只在「有東西壞掉時」增加可見性——但在大多
   數正常的週次，表格只會是 ~52 列幾乎相同的「✅ within ±0.05」。
3. 成本：新增一支腳本、一個工作流程工作、持續維護，以及在排程
   上增加 `contents: write` 權限。

**重新審視時機：** （a）已設定 `ALPHAGENOME_API_KEY` GitHub
   secret，且（b）我們希望有一個使用者可見的「活躍 vs 棄用」
   訊號。完整的試探設計已保存於
   `.considerations/live-atlas-docs-page.md`（位於 repo 根目錄；
   或搜尋 git 紀錄中的 v0.24.0 試探 commit）。

### 更多條目將陸續新增

凡是透過此模式拒絕的功能，請在關閉該試探的同一個 commit 中，
將裁決記錄到這裡。未來的讀者應該可以 grep CONTRIBUTING 中的
「考慮過但未建構」關鍵字，看到所有曾被試探但未出貨的功能請求
完整歷史。

## 授權

貢獻程式碼即代表您同意自己的貢獻將依本專案的雙重授權條款授權
（AGPL-3.0-or-later 用於開源使用；商業授權可來信索取）。
