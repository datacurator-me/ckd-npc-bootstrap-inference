# Automatically step up to the project root if run from inside 'scripts/' or 'Scripts/'
if (tolower(basename(getwd())) == "scripts") {
  setwd("..")
}

dir.create("dataset", showWarnings = FALSE)
csv_path <- file.path("dataset", "kidney_disease.csv")

if (!file.exists(csv_path)) {
  url <- "https://raw.githubusercontent.com/vamshikrishnajr/chronic-kidney-disease-diagnosis/master/kidney_disease.csv"
  download.file(url, destfile = csv_path, quiet = TRUE)
}

df <- read.csv(csv_path, stringsAsFactors = FALSE, na.strings = c("?", "\t?", "NA", ""))
df$classification <- trimws(df$classification)
df$classification[df$classification == "ckd\t"] <- "ckd"
df <- df[df$classification %in% c("ckd", "notckd"), ]
df$classification <- factor(df$classification, levels = c("notckd", "ckd"))

biomarkers <- c("sc", "bu", "hemo", "sg", "bp")
for (v in biomarkers) df[[v]] <- as.numeric(df[[v]])

impute_conditional_mean <- function(x, group) {
  means <- tapply(x, group, mean, na.rm = TRUE)
  ifelse(is.na(x), means[group], x)
}

for (v in biomarkers) {
  df[[v]] <- impute_conditional_mean(df[[v]], df$classification)
}

library(coin)
library(boot)

cat("\n--- REAL KAGGLE CKD DATASET (n = 400) ---\n")
for (v in biomarkers) {
  m_notckd <- mean(df[[v]][df$classification == "notckd"])
  s_notckd <- sd(df[[v]][df$classification == "notckd"])
  m_ckd <- mean(df[[v]][df$classification == "ckd"])
  s_ckd <- sd(df[[v]][df$classification == "ckd"])
  cat(sprintf("%-6s: CKD = %7.3f +/- %6.3f | Non-CKD = %7.3f +/- %5.3f | Diff = %+7.3f\n",
              v, m_ckd, s_ckd, m_notckd, s_notckd, m_ckd - m_notckd))
}

set.seed(123)
cat("\n--- PERMUTATION TESTS (B = 2000) ---\n")
p_vals <- numeric(length(biomarkers))
names(p_vals) <- biomarkers

for (i in seq_along(biomarkers)) {
  v <- biomarkers[i]
  fmla <- as.formula(paste(v, "~ classification"))
  pt <- independence_test(fmla, data = df, distribution = approximate(nresample = 2000))
  pv <- as.numeric(pvalue(pt))
  p_vals[i] <- max(pv, 1 / 2001)
  cat(sprintf("%-6s: p-value = %0.4f (Z = %6.2f)\n", v, p_vals[i], as.numeric(statistic(pt))))
}

lambda_obs <- -2 * sum(log(p_vals))
k <- length(biomarkers)
p_global <- 1 - pchisq(lambda_obs, df = 2 * k)
cat(sprintf("\nFisher Omnibus Lambda = %0.2f, df = %d, p_global = %0.4e\n", lambda_obs, 2 * k, p_global))

cat("\n--- BOOTSTRAP 95% CIs (B = 1000) ---\n")
mean_diff_fn <- function(data, idx, var_name) {
  d <- data[idx, ]
  mean(d[[var_name]][d$classification == "ckd"]) - mean(d[[var_name]][d$classification == "notckd"])
}

for (v in biomarkers) {
  b <- boot(data = df, statistic = function(d, idx) mean_diff_fn(d, idx, v), R = 1000)
  ci <- boot.ci(b, type = "perc", conf = 0.95)
  cat(sprintf("%-6s: Diff = %+7.3f, 95%% CI = [%+7.3f, %+7.3f], SE = %0.3f\n",
              v, b$t0, ci$percent[4], ci$percent[5], sd(b$t)))
}