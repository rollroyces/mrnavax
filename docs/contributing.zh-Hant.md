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
