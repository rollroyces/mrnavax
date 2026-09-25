# Atlas API 金鑰啟用

工具包預設帶有真實資料的 AlphaGenome Atlas 整合，但在即時每週排程可以拉取真實資料之前，需要先設定 Google AI Studio API 金鑰。本頁涵蓋**三種啟用路徑**：

## 為什麼有三種路徑？

金鑰最多可以存在三個位置：

  1. `ALPHAGENOME_API_KEY` **環境變數** — CI 偏好
  2. **`~/projects/alphagenome-work/.alphagenome_key`** — 本地開發偏好（自動偵測；v0.25.1 新增）
  3. **`~/.alphagenome_key`** — 替代 home 位置（v0.25.1 新增）

載入器依此順序解析：明確的 `api_key=...` 參數 > `ALPHAGENOME_API_KEY` 環境變數 > helper file `~/projects/alphagenome-work/.alphagenome_key` > `~/.alphagenome_key`。

## 路徑 1：本地開發（推薦）

如果你 clone 了這個 repo 並希望 Atlas 呼叫從筆電「直接可用」，請將 39 字元的 Google AI Studio 金鑰（以 `AIza` 開頭）寫入：

```bash
mkdir -p ~/projects/alphagenome-work
echo -n 'YOUR_39_CHAR_KEY_HERE' > ~/projects/alphagenome-work/.alphagenome_key
chmod 600 ~/projects/alphagenome-work/.alphagenome_key
```

下次執行 `score_variant()` 或 `live_atlas_integration.py` 時，載入器會自動偵測此檔案。

## 路徑 2：GitHub Actions 每週排程（CI 推薦）

要啟用每週的 Atlas 整合測試：

  1. 從 https://aistudio.google.com/apikey 取得 Google AI Studio API 金鑰（提供免費方案）。
  2. 開啟 https://github.com/rollroyces/mrnavax/settings/secrets/actions/new
  3. 名稱：`ALPHAGENOME_API_KEY`
  4. 值：貼上 39 字元的金鑰
  5. 點擊「Add secret」

下次排定的 cron（每週一 06:00 UTC）或手動 `gh workflow run atlas_integration.yml` 將會呼叫即時 Atlas API，並將結果與 `tests/fixtures/alphagenome_atlas_live.json` 中的 fixture 比對。如果上游 API 契約變更，整合測試會在一週內捕捉到。

## 路徑 3：CI / 腳本（環境變數）

對於暫時性環境（Docker、沒有 helper file 的 CI runner）：

```bash
export ALPHAGENOME_API_KEY='YOUR_39_CHAR_KEY_HERE'
```

## 驗證啟用

設定任何路徑後，執行驗證腳本以確認端到端 Atlas 連線：

```bash
python scripts/check_atlas_activation.py
```

預期輸出（hg38 的 BRAF V600E）：

```
[helper-file] OK: /Users/hermes/projects/alphagenome-work/.alphagenome_key (40 bytes, prefix matches Google AI key: True)
[package] alphagenome installed
[key] OK: resolved 39-char key (prefix matches Google AI: True)
[live-call] score=1.000, classification='high', is_coding=True, n_scorers=19
[fixture] recorded=1.000/'high', live=1.000/'high', matches=True
[result] PASS: live Atlas matches recorded fixture
[result] Activation verified — weekly Atlas cron will work
```

如果即時呼叫失敗，腳本會印出上游錯誤訊息並以非零狀態退出。

## 「已啟用」的定義

Atlas API 金鑰「已啟用」代表**全部三項**皆為真：

  1. ✅ 金鑰能透過上述三種路徑之一成功解析
  2. ✅ alphagenome 套件已安裝（`pip install 'mrnavax[variant-alphagenome]'`）
  3. ✅ 即時 BRAF V600E 呼叫回傳 `score≈1.0`、`classification='high'`、`is_coding=True`（與記錄的 fixture 吻合）

當三項皆為真時，GitHub Actions 每週排程將執行即時 Atlas 呼叫，而非以「ALPHAGENOME_API_KEY not configured」通知跳過。

## 隱私

API 金鑰絕不會被記錄、回顯，或由工具包寫入磁碟。驗證腳本僅列印其長度與前綴匹配狀態，而非金鑰本身。