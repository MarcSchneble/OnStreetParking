# set current path to this file
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# clear the workspace and restart session
rm(list = ls())

# sets local system language to English
Sys.setlocale("LC_ALL","English")

source("Functions.R")

library(dplyr)
library(lubridate)
library(ggplot2)
library(survival)
library(grDevices)
library(pROC)
library(pracma)
library(readr)
library(spatstat)
library(readxl)
library(scales)
library(frailtypack)
library(mgcv)
library(PRROC)

# read parking data ----

G <- readRDS("Data/network.rds")
intens.G <- density.lpp(unmark(G), sigma = 50)

data <- read_rds("Data/data_2019_clean.rds") %>%
  filter(StreetMarker %in% levels(G$data$marks),
         same.day == 1, irregular == 0, DurationMinutes2 > 0,
         h.start >= 8, h.start < 20,
         d.start >= 91, d.start <= 191) %>%
  mutate(SideOfStreetCode = factor(SideOfStreetCode, levels = unique(SideOfStreetCode)),
         StreetMarker = factor(StreetMarker, levels = levels(G$data$marks)),
         State = State - 1)

G$data <- G$data[G$data$marks %in% data$StreetMarker]
G$data$marks <- factor(G$data$marks, levels = unique(G$data$marks))
data$StreetMarker <- factor(data$StreetMarker, levels = levels(G$data$marks))

# distance matrix
N <- length(levels(data$StreetMarker))
distance <- matrix(NA, N, N)
for (i in 1:N) {
  G.i <- G
  G.i$data <- G$data[i, ]
  fundist <- distfun(G.i)
  distance[i, ] <- fundist(G)
}

#source("occupancy.R")
data <- read_rds("Data/data_frac.rds")
data$frac0[which(is.na(data$frac0))] <- 0
data$frac1[which(is.na(data$frac1))] <- 0

# model Lonsdale east ----
data.lonsdale.east <- filter(data, StreetName == "LONSDALE STREET",
                             BetweenStreet1 == "RUSSELL STREET" | BetweenStreet1 == "EXHIBITION STREET",
                             BetweenStreet2 == "EXHIBITION STREET" | BetweenStreet2 == "SPRING STREET") %>%
  mutate(SideOfStreetCode = factor(SideOfStreetCode, levels = unique(SideOfStreetCode)))

sm <- smooth.construct.ps.smooth.spec(s(m.start), data = data.lonsdale.east, knots = NULL)
X <- sm$X
X <- sweep(X, 2, colMeans(X))[, -1]
Z <- as.data.frame(X)
colnames(Z) <- paste0("m", 1:ncol(X))

data.lonsdale.east <- bind_cols(data.lonsdale.east, Z)
ind.m <- match(unique(data.lonsdale.east$m.start), data.lonsdale.east$m.start)
x <- data.lonsdale.east$m.start[ind.m]

model.weibull.0 <- survreg(Surv(pmin(DurationMinutes2, 60), 1*(DurationMinutes2 <= 60)) ~ weekday + frailty(StreetMarker) +
                             SideOfStreetCode + frac1 +
                             m1 + m2 + m3 + m4 + m5 + m6 + m7 + m8 + m9,
                           dist = "weibull", data = data.lonsdale.east %>% filter(State == 0), score = TRUE)

par.m <- as.numeric(tail(model.weibull.0$coefficients, ncol(X)))
x <- data.lonsdale.east$m.start[ind.m]
y <- as.vector(X[ind.m, ]%*%par.m)
limits <- smoothConfidence(par.m,
                           model.weibull.0$var[11:19, 11:19],
                           X[ind.m, ])
df <- tibble(x = x, y = -y/model.weibull.0$scale, y.lower = -limits$lower/model.weibull.0$scale,
             y.upper = -limits$upper/model.weibull.0$scale)
g0 <- ggplot(df) +
  geom_line(aes(x = x, y = y), col = "red") +
  geom_ribbon(aes(x = x, ymin = y.lower, ymax = y.upper), alpha = 0.5) +
  theme_bw() +
  scale_x_continuous(breaks = seq(360, 1200, 120), labels = c("6am", "8am", "10am", "12pm", "2pm", "4pm", "6pm", "8pm")) +
  labs(x = expression(paste(hour[t])), y = expression(paste(g[4*","*0](hour[t]))), parse = TRUE) +
  scale_y_continuous(limits = c(-0.8, 1))
pdf(file = "Plots/hour0.pdf", width = 4.5, height = 3)
print(g0)
dev.off()


model.weibull.1 <- survreg(Surv(pmin(DurationMinutes2, 60), 1*(DurationMinutes2 <= 60)) ~ weekday + frailty(StreetMarker) +
                             SideOfStreetCode + frac0 +
                             m1 + m2 + m3 + m4 + m5 + m6 + m7 + m8 + m9,
                           dist = "weibull", data = data.lonsdale.east %>% filter(State == 1), score = TRUE)

par.m <- as.numeric(tail(model.weibull.1$coefficients, ncol(X)))
y <- as.vector(X[ind.m, ]%*%par.m)
limits <- smoothConfidence(par.m,
                           model.weibull.1$var[11:19, 11:19],
                           X[ind.m, ])
df <- tibble(x = x, y = -y/model.weibull.1$scale, y.lower = -limits$lower/model.weibull.1$scale,
             y.upper = -limits$upper/model.weibull.1$scale)
g1 <- ggplot(df) +
  geom_line(aes(x = x, y = y), col = "red") +
  geom_ribbon(aes(x = x, ymin = y.lower, ymax = y.upper), alpha = 0.5) +
  theme_bw() +
  scale_x_continuous(breaks = seq(360, 1200, 120), labels = c("6am", "8am", "10am", "12pm", "2pm", "4pm", "6pm", "8pm")) +
  labs(x = expression(paste(hour[t])), y = expression(paste(g[4*","*1](hour[t]))), parse = TRUE) +
  scale_y_continuous(limits = c(-0.8, 1))
pdf(file = "Plots/hour1.pdf", width = 4.5, height = 3)
print(g1)
dev.off()



# prediction ----
R <- 100
dur <- 10
time <- as.POSIXct("2019-06-01 10:00:00", tz = "Australia/Melbourne")
observed <- prediction.det <- prediction.exp <- prediction.weibull <- prediction.weibull.wostar <- NULL
for (r in 1:R) {
  time.r <- time + days(sample(0:29, 1)) + seconds(runif(1, 0, 60*120))
  # simulate data point and compute distance function
  P <- rlpp(1, intens.G)
  fundist <- distfun(P)
  dist <- fundist(G)
  # data for distance
  ind.distance <- which(dist <= 300)
  markers.distance <- levels(data$StreetMarker)[ind.distance]
  ind.fit <- which(dist <= 250)
  if (length(ind.fit) > 10)  {
    print(r)
    markers.fit <- levels(data$StreetMarker)[ind.fit]
    pre <- get_occupancy2(time.r, data, ind.distance, duration = TRUE)
    post <- get_occupancy2(time.r + minutes(dur), data, ind.fit, duration = FALSE)
    # fit models
    data.r <- data %>% filter(StreetMarker %in% markers.fit,
                              difftime(time.r, DepartureTime, units = "days") > 0,
                              difftime(time.r, DepartureTime, units = "days") < 30) %>%
      mutate(SideOfStreetCode = factor(SideOfStreetCode, levels = unique(SideOfStreetCode)),
             StreetMarker = factor(StreetMarker, levels = unique(StreetMarker)))
    data0 <- data.r %>% filter(State == 0)
    data1 <- data.r %>% filter(State == 1)
    markers <- intersect(unique(data0$StreetMarker), unique(data1$StreetMarker))
    data0 <- data0 %>% filter(StreetMarker %in% markers) %>%
      mutate(SideOfStreetCode = factor(SideOfStreetCode, levels = unique(SideOfStreetCode)),
             StreetMarker = factor(StreetMarker, levels = unique(StreetMarker)))
    data1 <- data1 %>% filter(StreetMarker %in% markers) %>%
      mutate(SideOfStreetCode = factor(SideOfStreetCode, levels = unique(SideOfStreetCode)),
             StreetMarker = factor(StreetMarker, levels = unique(StreetMarker)))
    fmla0 <- Surv(pmin(DurationMinutes2, 60), 1*(DurationMinutes2 <= 60)) ~ weekday + frailty(StreetMarker) +
      factor(h.start) + frac1
    if (length(unique(data0$SideOfStreetCode)) > 1) {
      fmla0 <- update(fmla0, ~ .+ SideOfStreetCode)
    }
    fmla1 <- Surv(pmin(DurationMinutes2, 60), 1*(DurationMinutes2 <= 60)) ~ weekday + frailty(StreetMarker) +
      factor(h.start) + frac1
    if (length(unique(data1$SideOfStreetCode)) > 1) {
      fmla1 <- update(fmla1, ~ .+ SideOfStreetCode)
    }
    model.exp.0 <- survreg(fmla0, dist = "exponential", data = data0, score = TRUE)
    model.exp.1 <- survreg(fmla1, dist = "exponential", data = data1, score = TRUE)
    model.weibull.0 <- survreg(fmla0, dist = "weibull", data = data0, score = TRUE)
    model.weibull.1 <- survreg(fmla1, dist = "weibull", data = data1, score = TRUE)
    # get parameters from the model
    par.exp.0 <- get_par2(time.r, model.exp.0, data0, pre, 0, distance)
    par.exp.1 <- get_par2(time.r, model.exp.1, data1, pre, 1, distance)
    par.weibull.0 <- get_par2(time.r, model.weibull.0, data0, pre, 0, distance)
    par.weibull.1 <- get_par2(time.r, model.weibull.1, data1, pre, 1, distance)
    # predict occupancy
    occupancy.pre <- filter(pre, marker %in% markers)$occupancy
    occupancy.post <- filter(post, marker %in% markers)$occupancy
    d0 <- filter(pre, marker %in% markers)$duration
    ind.prediction <- which(is.element(occupancy.pre, c(0, 1)) & is.element(occupancy.post, c(0, 1)))
    prediction.exp <- c(prediction.exp, get_prediction_exp(occupancy.pre, par.exp.0, par.exp.1, dur)[ind.prediction])
    prediction.weibull <- c(prediction.weibull, get_prediction_weibull(occupancy.pre, d0, par.weibull.0, par.weibull.1, dur)[ind.prediction])
    d0 <- rep(0, length(d0))
    prediction.weibull.wostar <- c(prediction.weibull.wostar, get_prediction_weibull(occupancy.pre, d0, par.weibull.0, par.weibull.1, dur)[ind.prediction])
    prediction.det <- c(prediction.det, 1 - occupancy.pre[ind.prediction])
    observed <- c(observed, occupancy.post[ind.prediction])
  }
}

roc.exp <- pROC::roc(observed, prediction.exp, levels = c(1, 0))
roc.weibull <- pROC::roc(observed, prediction.weibull, levels = c(1, 0))
roc.weibull.wostar <- pROC::roc(observed, prediction.weibull.wostar, levels = c(1, 0))
roc.check <- pROC::roc(observed, prediction.det, levels = c(1, 0))

df.exp <- tibble(x = roc.exp$specificities, y = roc.exp$sensitivities, kind = "exp")
df.weibull <- tibble(x = roc.weibull$specificities, y = roc.weibull$sensitivities, kind = "weibull")
df.weibull.wostar <- tibble(x = roc.weibull.wostar$specificities,
                            y = roc.weibull.wostar$sensitivities, kind = "weibullwostar")
df.check <- tibble(x = roc.check$specificities, y = roc.check$sensitivities, kind = "check")
df.random <- tibble(x = c(0, 1), y = c(1, 0), kind = "random")
df <- bind_rows(df.weibull, df.weibull.wostar, df.exp) %>%
  mutate(kind = factor(kind, levels = c("weibull", "weibullwostar", "exp")))

g <- ggplot() +
  geom_ribbon(data = df %>% filter(kind == "weibull"), aes(x = 1-x, ymin = 0, ymax = y),
              alpha = 0.2) +
  geom_line(data = df, aes(x = 1-x, y = y, color = kind, linetype = kind)) +
  scale_color_hue(name = "Predictor",
                  labels = c(paste0("Semi-Markov\n(AUC = ", round(roc.weibull$auc, 3), ")"),
                             paste0("Markov\n(AUC = ", round(roc.exp$auc, 3), ")"),
                             paste0("Random\n(AUC = 0.5)"))) +
  labs(x = "Specifciity", y = "Senisitivity", linetype = "Predictor") +
  theme_bw() +
  theme(legend.position = "bottom", legend.key.width = unit(0.8, "cm")) +
  scale_x_continuous(breaks = seq(0, 1, 0.2), labels = seq(1, 0, -0.2)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  scale_linetype_manual(labels = c(paste0("Semi-Markov\n(AUC = ", round(roc.weibull$auc, 3), ")"),
                                   paste0("Markov\n(AUC = ", round(roc.exp$auc, 3), ")"),
                                   paste0("Random\n(AUC = 0.5)")),
                        values = c(1, 3, 5))

pdf(file = "Plots/roc_30min_afternoon.pdf", width = 5.5, height = 5.5)
print(g)
dev.off()



# modeling with different prediction horizons ----
dur <- seq(5, 60, 5)
AUC <- matrix(0, length(dur), 3)
R <- 100
daytime <- 10
for (d in 1:length(dur)) {
  time <- as.POSIXct("2019-06-01", tz = "Australia/Melbourne") + hours(daytime)
  observed <- prediction.det <- prediction.exp <- prediction.weibull <- prediction.weibull.wostar <- NULL
  for (r in 1:R) {
    time.r <- time + days(sample(0:29, 1)) + seconds(runif(1, 0, 60*120))
    # simulate data point and compute distance function
    P <- rlpp(1, intens.G)
    fundist <- distfun(P)
    dist <- fundist(G)
    # data for distance
    ind.distance <- which(dist <= 300)
    markers.distance <- levels(data$StreetMarker)[ind.distance]
    ind.fit <- which(dist <= 250)
    if (length(ind.fit) > 10)  {
      print(r)
      markers.fit <- levels(data$StreetMarker)[ind.fit]
      pre <- get_occupancy2(time.r, data, ind.distance, duration = TRUE)
      post <- get_occupancy2(time.r + minutes(dur[d]), data, ind.fit, duration = FALSE)
      # fit models
      data.r <- data %>% filter(StreetMarker %in% markers.fit,
                                difftime(time.r, DepartureTime, units = "days") > 0,
                                difftime(time.r, DepartureTime, units = "days") < 30) %>%
        mutate(SideOfStreetCode = factor(SideOfStreetCode, levels = unique(SideOfStreetCode)),
               StreetMarker = factor(StreetMarker, levels = unique(StreetMarker)))
      data0 <- data.r %>% filter(State == 0)
      data1 <- data.r %>% filter(State == 1)
      markers <- intersect(unique(data0$StreetMarker), unique(data1$StreetMarker))
      data0 <- data0 %>% filter(StreetMarker %in% markers) %>%
        mutate(SideOfStreetCode = factor(SideOfStreetCode, levels = unique(SideOfStreetCode)),
               StreetMarker = factor(StreetMarker, levels = unique(StreetMarker)))
      data1 <- data1 %>% filter(StreetMarker %in% markers) %>%
        mutate(SideOfStreetCode = factor(SideOfStreetCode, levels = unique(SideOfStreetCode)),
               StreetMarker = factor(StreetMarker, levels = unique(StreetMarker)))
      fmla0 <- Surv(pmin(DurationMinutes2, 60), 1*(DurationMinutes2 <= 60)) ~ weekday + frailty(StreetMarker) +
        factor(h.start) + frac1
      if (length(unique(data0$SideOfStreetCode)) > 1) {
        fmla0 <- update(fmla0, ~ .+ SideOfStreetCode)
      }
      fmla1 <- Surv(pmin(DurationMinutes2, 60), 1*(DurationMinutes2 <= 60)) ~ weekday + frailty(StreetMarker) +
        factor(h.start) + frac1
      if (length(unique(data1$SideOfStreetCode)) > 1) {
        fmla1 <- update(fmla1, ~ .+ SideOfStreetCode)
      }
      model.exp.0 <- survreg(fmla0, dist = "exponential", data = data0, score = TRUE)
      model.exp.1 <- survreg(fmla1, dist = "exponential", data = data1, score = TRUE)
      model.weibull.0 <- survreg(fmla0, dist = "weibull", data = data0, score = TRUE)
      model.weibull.1 <- survreg(fmla1, dist = "weibull", data = data1, score = TRUE)
      # get parameters from the model
      par.exp.0 <- get_par2(time.r, model.exp.0, data0, pre, 0, distance)
      par.exp.1 <- get_par2(time.r, model.exp.1, data1, pre, 1, distance)
      par.weibull.0 <- get_par2(time.r, model.weibull.0, data0, pre, 0, distance)
      par.weibull.1 <- get_par2(time.r, model.weibull.1, data1, pre, 1, distance)
      # predict occupancy
      occupancy.pre <- filter(pre, marker %in% markers)$occupancy
      occupancy.post <- filter(post, marker %in% markers)$occupancy
      d0 <- filter(pre, marker %in% markers)$duration
      ind.prediction <- which(is.element(occupancy.pre, c(0, 1)) & is.element(occupancy.post, c(0, 1)))
      prediction.exp <- c(prediction.exp,
                          get_prediction_exp(occupancy.pre, par.exp.0, par.exp.1, dur[d])[ind.prediction])
      prediction.weibull <- c(prediction.weibull,
                              get_prediction_weibull(occupancy.pre, d0, par.weibull.0, par.weibull.1, dur[d])[ind.prediction])
      prediction.weibull.wostar <- c(prediction.weibull.wostar,
                                     get_prediction_weibull(occupancy.pre, d0, par.weibull.0, par.weibull.1, dur[d], star = FALSE)[ind.prediction])
      prediction.det <- c(prediction.det, 1 - occupancy.pre[ind.prediction])
      observed <- c(observed, occupancy.post[ind.prediction])
    }
  }

  roc.exp <- pROC::roc(observed, prediction.exp, levels = c(1, 0))
  roc.weibull <- pROC::roc(observed, prediction.weibull, levels = c(1, 0))
  roc.weibull.wostar <- pROC::roc(observed, prediction.weibull.wostar, levels = c(1, 0))
  roc.check <- pROC::roc(observed, prediction.det, levels = c(1, 0))

  #AUC[d, ] <- c(roc.weibull$auc, roc.weibull.wostar$auc, roc.exp$auc)
  #saveRDS(AUC, file = "AUC.rds")

  df.exp <- tibble(x = roc.exp$specificities, y = roc.exp$sensitivities, kind = "exp")
  df.weibull <- tibble(x = roc.weibull$specificities, y = roc.weibull$sensitivities, kind = "weibull")
  df.weibull.wostar <- tibble(x = roc.weibull.wostar$specificities,
                              y = roc.weibull.wostar$sensitivities, kind = "weibullwostar")
  df.check <- tibble(x = roc.check$specificities, y = roc.check$sensitivities, kind = "check")
  df.random <- tibble(x = c(0, 1), y = c(1, 0), kind = "random")
  df <- bind_rows(df.weibull, df.weibull.wostar, df.exp) %>%
    mutate(kind = factor(kind, levels = c("weibull", "weibullwostar", "exp")))

  g <- ggplot() +
    geom_ribbon(data = df %>% filter(kind == "weibull"), aes(x = 1-x, ymin = 0, ymax = y),
                alpha = 0.2) +
    geom_line(data = df, aes(x = 1-x, y = y, color = kind, linetype = kind)) +
    scale_color_hue(name = "Predictor",
                    labels = c(paste0("Semi-Markov\n(state space S*,\nAUC = ", round(roc.weibull$auc, 3), ")"),
                               paste0("Semi-Markov\n(state space S,\nAUC = ", round(roc.weibull.wostar$auc, 3), ")"),
                               paste0("Markov\n(state space S,\nAUC = ", round(roc.exp$auc, 3), ")"))) +
    labs(x = "Specifciity", y = "Senisitivity", linetype = "Predictor") +
    theme_bw() +
    theme(legend.position = "bottom", legend.key.width = unit(0.8, "cm")) +
    scale_x_continuous(breaks = seq(0, 1, 0.2), labels = seq(1, 0, -0.2)) +
    scale_y_continuous(breaks = seq(0, 1, 0.2)) +
    scale_linetype_manual(labels = c(paste0("Semi-Markov\n(state space S*,\nAUC = ", round(roc.weibull$auc, 3), ")"),
                                     paste0("Semi-Markov\n(state space S,\nAUC = ", round(roc.weibull.wostar$auc, 3), ")"),
                                     paste0("Markov\n(state space S,\nAUC = ", round(roc.exp$auc, 3), ")")),
                          values = c(1, 3, 5)) +
    geom_line(data = tibble(x = c(0, 1), y = c(0, 1)), aes(x = x, y = y))

  pdf(file = paste0("Plots/roc_", dur[d], "min_daytime", daytime, ".pdf"), width = 5.5, height = 5.5)
  print(g)
  dev.off()
}
AUC <- readRDS("AUC.rds")

df <- tibble(duration = rep(dur, 3),
             AUC = as.vector(AUC),
             predictor = rep(c("weibull", "weibullwostar", "exp"), each = length(dur))) %>%
  mutate(predictor = factor(predictor, levels = c("weibull", "weibullwostar", "exp")))
g <- ggplot(df, aes(x = duration, y = AUC, color = predictor)) +
  geom_point() +
  geom_line(aes(linetype = predictor)) +
  theme_bw() +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_x_continuous(breaks = seq(0, 60, 10)) +
  labs(x = "Prediction horizon (in minutes)", color = "Predictor", linetype = "Predictor") +
  scale_linetype_manual(labels = c(paste0("Semi-Markov (state space S*)"),
                                   paste0("Semi-Markov (state space S)"),
                                   paste0("Markov (State space S)")),
                        values = c(1, 3, 5)) +
  scale_color_hue(labels = c(paste0("Semi-Markov (state space S*)"),
                                      paste0("Semi-Markov (state space S)"),
                                      paste0("Markov (State space S)"))) +
  geom_hline(yintercept = 0.5)

pdf(file = paste0("Plots/AUC.pdf"), width = 8, height = 4)
print(g)
dev.off()

