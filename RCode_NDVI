## ====================================================================================================
## Post-fire NDVI recovery, Jonkershoek (fire of 21 Feb 2021)
## Exponential recovery model without a seasonal term:  NDVI(t) = alpha + gamma * (1 - e^(-t / lambda))
## ====================================================================================================

## ---- 0. Setup -----------------------------------------------------------------
## Requires: tidyverse, patchwork, and MASS (MASS comes with R)
## install.packages(c("tidyverse", "patchwork"))
library(tidyverse)


## ---- 1. Load data -------------------------------------------------------------
## Semicolon-separated, decimal POINTS (only the unused time_start column has commas). Check str(ndvi): NDVI and years_since_fire must be "num".
ndvi <- read.csv("Jonkershoek_NDVI_quarterly.csv", sep = ";",
                 stringsAsFactors = FALSE)
ndvi$date <- as.Date(ndvi$date, format = "%Y/%m/%d")
str(ndvi)
names(ndvi)

## >>> EDIT these two to match the column names printed by names(ndvi) <<<
ndvi_col <- "NDVI"
t_col    <- "years_since_fire"

stopifnot(is.numeric(ndvi[[ndvi_col]]), is.numeric(ndvi[[t_col]]))

pre  <- ndvi[ndvi$period == "pre-fire",  ]
post <- ndvi[ndvi$period == "post-fire", ]

dat <- data.frame(t = post[[t_col]], y = post[[ndvi_col]], date = post$date)
dat <- dat[complete.cases(dat), ]
dat <- dat[order(dat$t), ]
nrow(dat)                      # number of post-fire observations
summary(dat)

pre_ndvi <- mean(pre[[ndvi_col]], na.rm = TRUE)   # pre-fire baseline
pre_ndvi

## ---- 2. Quick look at the data ------------------------------------------------
ggplot(dat, aes(t, y)) +
  geom_point() +
  geom_hline(yintercept = pre_ndvi, linetype = "dashed") +
  labs(x = "Time since fire (years)", y = "Mean NDVI")

# day of year (1-366) from the Date column
dat$doy <- as.numeric(format(dat$date, "%j"))

# angle around the annual cycle
dat$w <- 2 * pi * dat$doy / 365.25

head(dat)
range(dat$doy)


## ---- 3. The linear model (model 1) --------------------------------------------
fit_lin <- lm(y ~ t, data = dat)

summary(fit_lin)
AIC(fit_lin)

ggplot(dat, aes(t, y)) +
  geom_point() +
  geom_line(aes(y = fitted(fit_lin)), colour = "red")

## ---- 4. Exponential recovery, no seasonality (model 2) ------------------------
start_exp <- list(alpha = min(dat$y),          # lowest observed y
                  gam   = max(dat$y) - min(dat$y),          # range of y: max minus min
                  lambda = 1)

fit_exp <- nls(y ~ alpha + gam * (1 - exp(-t / lambda)),               # write the model using alpha, gam, lambda, t
               data  = dat,
               start = start_exp)

summary(fit_exp)
AIC(fit_exp)
confint(fit_exp)       # 95% CIs for your results table (may need library(MASS) on older R)

## Add the fitted curve to the plot, next to the linear model
ggplot(dat, aes(t, y)) +
  geom_point() +
  geom_hline(yintercept = pre_ndvi, linetype = "dashed") +
  geom_line(aes(y = fitted(fit_lin)), colour = "red") +
  geom_line(aes(y = fitted(fit_exp)), colour = "blue") +
  labs(x = "Time since fire (years)", y = "Mean NDVI")

## ---- 5. Full model, exponential recovery + seasonality (model 3) --------------
## NDVI(t) = alpha + gam*(1 - exp(-t/lambda)) + b1*sin(w) + b2*cos(w)
dat$s <- sin(dat$w)
dat$c <- cos(dat$w)

start_sea <- list(alpha  = min(dat$y),
                  gam    = max(dat$y) - min(dat$y),
                  lambda = 1,
                  b1     = 0,
                  b2     = 0)

fit_sea <- nls(y ~ alpha + gam * (1 - exp(-t / lambda)) + b1 * s + b2 * c,
               data  = dat,
               start = start_sea)

summary(fit_sea)
confint(fit_sea)

## Recover amplitude and phase from b1, b2
b1 <- coef(fit_sea)["b1"]
b2 <- coef(fit_sea)["b2"]
A   <- sqrt(b1^2 + b2^2)
phi <- atan2(b2, b1)
A; phi

## Plateau and recovery time
alpha  <- coef(fit_sea)["alpha"]
gam    <- coef(fit_sea)["gam"]
lambda <- coef(fit_sea)["lambda"]
plateau <- alpha + gam
RT      <- lambda * log(200 * gam)
plateau; RT

## Compare all three models
AIC(fit_lin, fit_exp, fit_sea)

## Plot all three
ggplot(dat, aes(t, y)) +
  geom_point() +
  geom_hline(yintercept = pre_ndvi, linetype = "dashed") +
  geom_line(aes(y = fitted(fit_lin)), colour = "red") +
  geom_line(aes(y = fitted(fit_exp)), colour = "blue") +
  geom_line(aes(y = fitted(fit_sea)), colour = "darkgreen", linewidth = 1) +
  labs(x = "Time since fire (years)", y = "Mean NDVI")

## ---- 6.a Residual diagnostics for the full model--------------------------------
dat$res    <- resid(fit_sea)
dat$fitted <- fitted(fit_sea)

## (a) Residuals over time: any leftover pattern?
ggplot(dat, aes(t, res)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Time since fire (years)", y = "Residual (NDVI)")

## (b) Residuals vs fitted: is the spread constant?
ggplot(dat, aes(fitted, res)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Fitted NDVI", y = "Residual (NDVI)")

## (c) Normality
qqnorm(dat$res); qqline(dat$res)
shapiro.test(dat$res)

## (d) Autocorrelation: are neighbouring residuals related?
acf(dat$res, main = "ACF of residuals")

## Biggest misses
dat[order(-abs(dat$res)), ][1:5, c("t", "date", "y", "fitted", "res")]

## ---- 6.b All residual diagnostics on one plot----------------------------------
library(patchwork)

dat$res    <- resid(fit_sea)
dat$fitted <- fitted(fit_sea)
sw <- shapiro.test(dat$res)

## (a) Residuals over time
pa <- ggplot(dat, aes(t, res)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Time since fire (years)", y = "Residual (NDVI)",
       title = "Residuals over time")

## (b) Residuals vs fitted
pb <- ggplot(dat, aes(fitted, res)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Fitted NDVI", y = "Residual (NDVI)",
       title = "Residuals vs fitted")

## (c) Normal Q-Q
pc <- ggplot(dat, aes(sample = res)) +
  stat_qq() +
  stat_qq_line() +
  labs(x = "Theoretical quantiles", y = "Sample quantiles",
       title = "Normal Q-Q",
       subtitle = paste0("Shapiro-Wilk p = ", signif(sw$p.value, 2)))

## (d) Autocorrelation of residuals
a      <- acf(dat$res, plot = FALSE)
acf_df <- data.frame(lag = a$lag[-1], acf = a$acf[-1])
bound  <- qnorm(0.975) / sqrt(nrow(dat))

pd <- ggplot(acf_df, aes(lag, acf)) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = c(-bound, bound), linetype = "dashed", colour = "blue") +
  geom_segment(aes(xend = lag, yend = 0)) +
  labs(x = "Lag", y = "ACF", title = "Autocorrelation of residuals")

## Combine into a 2 x 2 grid with A-D labels
diag_fig <- (pa | pb) / (pc | pd) +
  plot_annotation(tag_levels = "A") &
  theme_bw(base_size = 11)

diag_fig
ggsave("Figure_S1_residual_diagnostics.png", diag_fig,
       width = 9, height = 7, dpi = 300)

## Residual-versus-time and residuals-versus-fitted plots
pa <- ggplot(dat, aes(t, res)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Time since fire (years)", y = "Residual (NDVI)")

pb <- ggplot(dat, aes(fitted, res)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Fitted NDVI", y = "Residual (NDVI)")

res_fig <- (pa | pb) + plot_annotation(tag_levels = "A") & theme_bw()
res_fig

ggsave("Figure_S1_residuals.png", res_fig, width = 9, height = 4, dpi = 300)

## ---- 7. Sensitvity to the 5 largest negative residuals ------------------------
drop_idx <- order(dat$res)[1:5]          # the 5 most negative residuals
dat_s    <- dat[-drop_idx, ]

fit_sea_s <- nls(y ~ alpha + gam * (1 - exp(-t / lambda)) + b1 * s + b2 * c,
                 data  = dat_s,
                 start = start_sea)

summary(fit_sea_s)
shapiro.test(resid(fit_sea_s))

## Compare key quantities, full data vs. 5 points removed
get_summary <- function(f) {
  p <- coef(f)
  c(alpha   = unname(p["alpha"]),
    gam     = unname(p["gam"]),
    lambda  = unname(p["lambda"]),
    A       = unname(sqrt(p["b1"]^2 + p["b2"]^2)),
    plateau = unname(p["alpha"] + p["gam"]),
    RT      = unname(p["lambda"] * log(200 * p["gam"])))
}
round(rbind(full = get_summary(fit_sea),
            drop5 = get_summary(fit_sea_s)), 4)

## ---- 8. Results tables for the report -----------------------------------------
## 8a. Parameter table (full model)
sm <- summary(fit_sea)$coefficients
ci <- confint(fit_sea)
tab_par <- data.frame(Parameter = rownames(sm),
                      Estimate  = sm[, "Estimate"],
                      SE        = sm[, "Std. Error"],
                      CI_low    = ci[, 1],
                      CI_high   = ci[, 2],
                      row.names = NULL)
tab_par[, -1] <- round(tab_par[, -1], 4)
tab_par

## 8b. Model comparison
aic_tab <- AIC(fit_lin, fit_exp, fit_sea)
aic_tab$dAIC   <- aic_tab$AIC - min(aic_tab$AIC)
aic_tab$weight <- exp(-aic_tab$dAIC / 2) / sum(exp(-aic_tab$dAIC / 2))
round(aic_tab, 3)

## 8c. Derived quantities with 95% intervals by simulation
##     (draw parameter sets from the fitted model's covariance matrix)
set.seed(1)
sims <- MASS::mvrnorm(10000, mu = coef(fit_sea), Sigma = vcov(fit_sea))
sims <- sims[sims[, "gam"] > 0.005 & sims[, "lambda"] > 0, ]   # RT needs gam > 0.005

derive <- function(p) {
  p   <- matrix(p, ncol = 5, dimnames = list(NULL, names(coef(fit_sea))))
  phi <- atan2(p[, "b2"], p[, "b1"])
  cbind(A         = sqrt(p[, "b1"]^2 + p[, "b2"]^2),
        plateau   = p[, "alpha"] + p[, "gam"],
        half_life = p[, "lambda"] * log(2),
        RT        = p[, "lambda"] * log(200 * p[, "gam"]),
        peak_doy  = ((pi / 2 - phi) %% (2 * pi)) / (2 * pi) * 365.25)
}

est  <- derive(coef(fit_sea))
sims_d <- derive(sims)
tab_der <- data.frame(Quantity = colnames(est),
                      Estimate = round(as.numeric(est), 3),
                      CI_low   = round(apply(sims_d, 2, quantile, 0.025), 3),
                      CI_high  = round(apply(sims_d, 2, quantile, 0.975), 3),
                      row.names = NULL)
tab_der

citation("tidyverse")
citation("patchwork")

