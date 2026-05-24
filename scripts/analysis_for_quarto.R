# scripts/analysis_for_quarto.R
# Quarto から source して実行するための単一ファイル（tidyverse を使用）
# 実行例: source("scripts/analysis_for_quarto.R")

library(tidyverse)
library(lubridate)

# 出力フォルダ（Quarto のルートからの相対パス）
out_dir <- "outputs"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# データ URL（TidyTuesday）
water_url   <- "https://raw.githubusercontent.com/rfordatascience/tidytuesday/main/data/2025/2025-05-20/water_quality.csv"
weather_url <- "https://raw.githubusercontent.com/rfordatascience/tidytuesday/main/data/2025/2025-05-20/weather.csv"

# ----- データ読み込み -----
cat("データを読み込みます...\n")
water <- readr::read_csv(water_url, show_col_types = FALSE) %>%
  mutate(date = ymd(date),
         enterococci_cfu_100ml = as.numeric(enterococci_cfu_100ml))
weather <- readr::read_csv(weather_url, show_col_types = FALSE) %>%
  mutate(date = ymd(date),
         precipitation_mm = as.numeric(precipitation_mm))

# ----- 天気を日次に集計 -----
weather_daily <- weather %>%
  group_by(date) %>%
  summarise(precipitation_mm = sum(precipitation_mm, na.rm = TRUE), .groups = "drop") %>%
  arrange(date)

# ----- ラグ付き降水（過去7日）を計算するヘルパー -----
compute_lagged <- function(dates, k, weather_df) {
  sapply(dates, function(d) {
    if (is.na(d)) return(NA_real_)
    start <- d - days(k)
    sum(weather_df$precipitation_mm[weather_df$date > start & weather_df$date <= d], na.rm = TRUE)
  })
}

unique_dates <- sort(unique(c(water$date, weather_daily$date)))
lagged_tbl <- tibble(date = unique_dates,
                     precip_7d = compute_lagged(unique_dates, 7, weather_daily))

# ----- データ結合 -----
data <- water %>%
  left_join(lagged_tbl, by = "date") %>%
  mutate(year = year(date))

# ----- 年次サマリー -----
yearly_summary <- data %>%
  group_by(year) %>%
  summarise(mean_enterococci = mean(enterococci_cfu_100ml, na.rm = TRUE),
            median_enterococci = median(enterococci_cfu_100ml, na.rm = TRUE),
            n = sum(!is.na(enterococci_cfu_100ml)),
            .groups = "drop") %>%
  arrange(year)

# 保存（CSV）
readr::write_csv(yearly_summary, file.path(out_dir, "yearly_summary.csv"))

# ----- 上位サイト（サンプル多いサイト）の年次推移（表示用） -----
top_sites <- data %>% count(swim_site, sort = TRUE) %>% slice_head(n = 6) %>% pull(swim_site)
site_yearly <- data %>%
  filter(swim_site %in% top_sites) %>%
  group_by(swim_site, year) %>%
  summarise(mean_enterococci = mean(enterococci_cfu_100ml, na.rm = TRUE), .groups = "drop")

# ----- プロットの準備（NA や非有限値を除去） -----
data_clean <- data %>%
  filter(!is.na(enterococci_cfu_100ml) & is.finite(enterococci_cfu_100ml) &
         !is.na(precip_7d) & is.finite(precip_7d))

library(ggplot2)

p1 <- ggplot(yearly_summary, aes(x = year, y = mean_enterococci)) +
  geom_line() + geom_point() +
  labs(title = "年ごとの平均 Enterococci", x = "年", y = "平均 Enterococci (CFU/100ml)") +
  theme_minimal()

p2 <- ggplot(site_yearly, aes(x = year, y = mean_enterococci, color = swim_site)) +
  geom_line() + geom_point() +
  labs(title = "主要泳場の年次平均 Enterococci", x = "年", y = "平均 Enterococci") +
  theme_minimal()

p3 <- ggplot(data_clean %>% mutate(rain_cat = ntile(precip_7d, 3)) %>% filter(!is.na(rain_cat)),
             aes(x = factor(rain_cat), y = enterococci_cfu_100ml)) +
  geom_boxplot() +
  labs(title = "降水（過去7日）カテゴリ別の Enterococci 分布",
       x = "降水カテゴリ（1=低,3=高）", y = "Enterococci (CFU/100ml)") +
  theme_minimal()

p4 <- ggplot(data_clean, aes(x = precip_7d, y = enterococci_cfu_100ml)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", se = FALSE) +   # 安定のため線形回帰を使用
  labs(title = "降水（過去7日）と Enterococci（散布図 + 回帰線）",
       x = "7日間累積降水 (mm)", y = "Enterococci (CFU/100ml)") +
  theme_minimal()

# ----- プロットを保存 -----
cat("プロットを保存します...\n")
ggsave(file.path(out_dir, "yearly_mean_enterococci.png"), p1, width = 8, height = 4)
ggsave(file.path(out_dir, "top_sites_yearly.png"), p2, width = 8, height = 5)
ggsave(file.path(out_dir, "enterococci_by_rain_cat.png"), p3, width = 6, height = 4)
ggsave(file.path(out_dir, "enterococci_vs_precip7d.png"), p4, width = 6, height = 4)

# ----- 図から読み取れるメッセージを自動生成 -----
# 1) 年次傾向（線形検定）
lm_year <- tryCatch(lm(mean_enterococci ~ year, data = yearly_summary), error = function(e) NULL)
msg_year <- if (is.null(lm_year)) {
  "年次傾向の検定を実行できませんでした。"
} else {
  s <- coef(summary(lm_year))["year", "Estimate"]
  p <- coef(summary(lm_year))["year", "Pr(>|t|)"]
  if (!is.na(p) && p < 0.05) {
    if (s > 0) {
      sprintf("年次平均 Enterococci は有意に増加しています（傾き=%.3f, p=%.3g）。", s, p)
    } else {
      sprintf("年次平均 Enterococci は有意に減少しています（傾き=%.3f, p=%.3g）。", s, p)
    }
  } else {
    sprintf("明確な線形トレンドは検出されません（p=%.3g）。年ごとの変動は存在します。", p)
  }
}

# 2) サイト別傾向（上位サイトの傾き上位を表示）
site_slopes <- data %>%
  filter(!is.na(enterococci_cfu_100ml)) %>%
  group_by(swim_site, year) %>%
  summarise(mean_e = mean(enterococci_cfu_100ml, na.rm = TRUE), .groups = "drop") %>%
  group_by(swim_site) %>%
  filter(n() >= 3) %>%
  summarise(slope = tryCatch(coef(lm(mean_e ~ year))["year"], error = function(e) NA_real_)) %>%
  arrange(desc(slope))

top_increasing_sites <- head(site_slopes$swim_site, 5)
msg_sites <- if (length(top_increasing_sites) == 0) {
  "十分なデータのあるサイトがありません。"
} else {
  paste0("増加傾向が大きいサイト例: ", paste(top_increasing_sites, collapse = ", "), "。個別調査を推奨します。")
}

# 3) 降水と Enterococci の相関（スピアマン）
full_clean <- data_clean
msg_precip <- "データ不足により相関検定を行えませんでした。"
if (nrow(full_clean) > 10) {
  ct <- tryCatch(cor.test(full_clean$precip_7d, full_clean$enterococci_cfu_100ml, method = "spearman"), error = function(e) NULL)
  if (!is.null(ct)) {
    rho <- ct$estimate; pval <- ct$p.value
    if (pval < 0.05 && abs(rho) > 0.1) {
      msg_precip <- sprintf("7日累積降水と Enterococci に有意な相関があります（rho=%.3f, p=%.3g）。降水増加後の濃度上昇が示唆されます（因果ではありません）。", rho, pval)
    } else {
      msg_precip <- sprintf("7日累積降水と Enterococci の相関は弱いか有意ではありません（rho=%.3f, p=%.3g）。", rho, pval)
    }
  }
}

# 結果をコンソール表示と CSV 保存
cat("------ 自動生成メッセージ ------\n")
cat(msg_year, "\n")
cat(msg_sites, "\n")
cat(msg_precip, "\n")

messages <- tibble(
  item = c("year_message", "site_message", "precip_message"),
  text = c(msg_year, msg_sites, msg_precip)
)
readr::write_csv(messages, file.path(out_dir, "report_messages.csv"))

cat("完了: outputs/ に図と report_messages.csv を保存しました。\n")