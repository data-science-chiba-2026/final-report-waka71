# Sydney Water Quality - 初心者向け解析スクリプト
# このスクリプトは tidyverse を使ってデータを読み込み、前処理し、
# 降水（過去1〜7日）の影響を調べ、グラフと概要CSVを outputs/ に保存します。

# 必要なパッケージ
# install.packages(c('tidyverse','lubridate'))  # 必要なら最初に実行してください
library(tidyverse)
library(lubridate)

# データのURL（TidyTuesday の生データ）
water_url <- "https://raw.githubusercontent.com/rfordatascience/tidytuesday/main/data/2025/2025-05-20/water_quality.csv"
weather_url <- "https://raw.githubusercontent.com/rfordatascience/tidytuesday/main/data/2025/2025-05-20/weather.csv"

# 出力フォルダ（scripts/から実行する想定）
out_dir <- file.path("..", "outputs")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# --- 1) データ読み込み ---
cat("データをダウンロードしています...\n")
water <- readr::read_csv(water_url, show_col_types = FALSE)
weather <- readr::read_csv(weather_url, show_col_types = FALSE)

# --- 2) 簡単な確認 ---
cat("水質データの先頭行:\n")
print(head(water))
cat("天気データの先頭行:\n")
print(head(weather))

# --- 3) 前処理 ---
# 日付を解析し、必要な列を数値化する
water <- water %>%
  mutate(
    date = lubridate::ymd(date),
    enterococci_cfu_100ml = as.numeric(enterococci_cfu_100ml),
    water_temperature_c = as.numeric(water_temperature_c),
    conductivity_ms_cm = as.numeric(conductivity_ms_cm)
  )

weather <- weather %>%
  mutate(
    date = lubridate::ymd(date),
    precipitation_mm = as.numeric(precipitation_mm)
  )

# 天気データを日付ごとに集計
weather_daily <- weather %>%
  group_by(date) %>%
  summarise(precipitation_mm = sum(precipitation_mm, na.rm = TRUE), .groups = "drop") %>%
  arrange(date)

# --- 4) ラグ付き降水量の計算（過去1〜7日） ---
compute_lagged <- function(dates, k, weather_df) {
  sapply(dates, function(d) {
    if (is.na(d)) return(NA_real_)
    start <- d - days(k)
    sum(weather_df$precipitation_mm[weather_df$date > start & weather_df$date <= d], na.rm = TRUE)
  })
}

unique_dates <- sort(unique(c(water$date, weather_daily$date)))
lagged_tbl <- tibble(date = unique_dates)
for (k in 1:7) {
  colname <- paste0("precip_lag", k, "d")
  lagged_tbl[[colname]] <- compute_lagged(lagged_tbl$date, k, weather_daily)
}
lagged_tbl <- lagged_tbl %>% rename(precip_7d = precip_lag7d)

# --- 5) データ結合 ---
data <- water %>% left_join(lagged_tbl %>% select(date, precip_7d), by = "date")

# --- 6) 集約と要約 ---
data <- data %>% mutate(year = lubridate::year(date))

yearly_summary <- data %>%
  group_by(year) %>%
  summarise(
    mean_enterococci = mean(enterococci_cfu_100ml, na.rm = TRUE),
    median_enterococci = median(enterococci_cfu_100ml, na.rm = TRUE),
    n = sum(!is.na(enterococci_cfu_100ml)),
    .groups = "drop"
  )

# 上位サイトの抽出
library(dplyr)
top_sites <- data %>% count(swim_site, sort = TRUE) %>% slice_head(n = 6) %>% pull(swim_site)
site_yearly <- data %>% filter(swim_site %in% top_sites) %>% group_by(swim_site, year) %>% summarise(mean_enterococci = mean(enterococci_cfu_100ml, na.rm = TRUE), .groups = "drop")

# 降水カテゴリ（3分位）
data <- data %>% mutate(rain_cat = ntile(precip_7d, 3)) %>% mutate(rain_cat = case_when(is.na(precip_7d) ~ "no_data", rain_cat == 1 ~ "low", rain_cat == 2 ~ "medium", rain_cat == 3 ~ "high"))

# プロット用に NA/無限大を除去したデータを作成
data_clean <- data %>% filter(!is.na(enterococci_cfu_100ml) & is.finite(enterococci_cfu_100ml) & !is.na(precip_7d) & is.finite(precip_7d))

# --- 7) プロット ---
library(ggplot2)

p1 <- ggplot(yearly_summary, aes(x = year, y = mean_enterococci)) + geom_line() + geom_point() + labs(title = "年ごとの平均 Enterococci", x = "年", y = "平均 Enterococci (CFU/100ml)") + theme_minimal()

p2 <- ggplot(site_yearly, aes(x = year, y = mean_enterococci, color = swim_site)) + geom_line() + geom_point() + labs(title = "主要泳場の年次平均 Enterococci", x = "年", y = "平均 Enterococci") + theme_minimal()

p3 <- ggplot(data %>% filter(!is.na(rain_cat) & rain_cat != "no_data"), aes(x = rain_cat, y = enterococci_cfu_100ml)) + geom_boxplot() + scale_y_continuous(trans = 'pseudo_log') + labs(title = "降水量（過去7日）カテゴリ別の Enterococci 分布", x = "降水カテゴリ", y = "Enterococci (CFU/100ml)") + theme_minimal()

p4 <- ggplot(data, aes(x = precip_7d, y = enterococci_cfu_100ml)) + geom_point(alpha = 0.3) + geom_smooth(method = "loess", se = TRUE) + labs(title = "降水（過去7日）と Enterococci", x = "7日間累積降水 (mm)", y = "Enterococci (CFU/100ml)") + theme_minimal()

# --- 8) 保存 ---
cat("プロットを保存しています...\n")
tryCatch({
  ggsave(filename = file.path(out_dir, "yearly_mean_enterococci.png"), plot = p1, width = 8, height = 4)
  ggsave(filename = file.path(out_dir, "top_sites_yearly.png"), plot = p2, width = 8, height = 5)
  ggsave(filename = file.path(out_dir, "enterococci_by_rain_cat.png"), plot = p3, width = 6, height = 4)
  ggsave(filename = file.path(out_dir, "enterococci_vs_precip7d.png"), plot = p4, width = 6, height = 4)
}, error = function(e) cat("プロット保存中にエラー: ", e$message, "\n"))

cat("要約CSVを保存しています...\n")
readr::write_csv(yearly_summary, file.path(out_dir, "yearly_summary.csv"))

cat("解析が完了しました。outputs/ フォルダに画像とCSVが保存されています。\n")
