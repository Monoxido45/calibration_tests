library(tidyverse)
# binary calibration package
library(CalibratR)
# graphical packages to reliability diagrams
library(grid)
library(ggpubr)

# reliability diagram, ECE and MCE
reliab_plot <- function(preds, resp, nbins = 50, title = ""){
  x_grid <- seq(0, 1, length.out = nbins)
  b_m <- x_grid[2] - x_grid[1]
  bin_data <- data.frame(pred_bin = preds,
                         resp_bin = resp)
  
  
  bins_stats <- x_grid[-1] |>
    map_dfr(~bin_data |>
              filter(pred_bin <= .x & pred_bin > (.x - b_m))|>
              summarise(acc = mean((resp_bin == 1)),
                        conf = mean(pred_bin)) |>
              mutate(gap = abs(acc - conf)) |>
              mutate(lim_sup =.x))
  
  bins_stats |>
    ggplot(aes(x = conf, y = acc)) +
    geom_point(colour = "dodgerblue3",
               alpha = 0.5) +
    geom_line(colour = "dodgerblue3") +
    geom_linerange(aes(x = conf, ymin = conf, 
                       ymax = acc), width = .01)+
    geom_abline(intercept = 0,
                linetype = "dashed",
                colour = "red") +
    theme_bw() +
    labs(y = "Calibration Function",
         x = "Probability estimate",
         title = title)+
    scale_x_continuous(breaks = scales::pretty_breaks(8),
                       limits = c(-0.05,1.05))+
    scale_y_continuous(breaks = scales::pretty_breaks(8),
                       limits = c(-0.05,1.05))
}

ece <- function(prob_estim, ground_truth, nbins = 50){
  x_grid <- seq(0, 1, length.out = nbins)
  b_m <- x_grid[2] - x_grid[1]
  bin_data <- data.frame(pred_bin = prob_estim,
                         resp_bin = ground_truth)
  P = length(ground_truth)
  
  x_grid[-1] |>
    map_dfr(~bin_data |>
              filter(pred_bin <= .x & pred_bin > (.x - b_m)) |>
              summarise(acc = mean((resp_bin == 1)),
                        conf = mean(pred_bin),
                        B = n()) |>
              mutate(gap = abs(acc - conf),
                     ECE = (B/P)*gap)) |>
    summarise(ECE = sum(ECE, na.rm = TRUE)) |>
    pull(ECE)
}

mce <- function(prob_estim, ground_truth, nbins = 50){
  x_grid <- seq(0, 1, length.out = nbins)
  b_m <- x_grid[2] - x_grid[1]
  bin_data <- data.frame(pred_bin = prob_estim,
                         resp_bin = ground_truth)
  P = length(ground_truth)
  
  x_grid[-1] |>
    map_dfr(~bin_data |>
              filter(pred_bin <= .x & pred_bin > (.x - b_m)) |>
              summarise(acc = mean((resp_bin == 1)),
                        conf = mean(pred_bin),
                        B = n()) |>
              mutate(gap = abs(acc - conf))) |>
    summarise(MCE = max(gap, na.rm = TRUE)) |>
    pull(MCE)
}

# Amazon reviews dataset for some examples
set.seed(1)
library(data.table)
library(tm)
library(glmnet)

dados <-  fread("data/Reviews.csv",
                header = TRUE,
                nrows = 5000) |>
  mutate(score_fat = ifelse(Score >= 4, 1, 0))

corp <- VCorpus(VectorSource(dados$Text))
dtm <- DocumentTermMatrix(corp,
                          control = list(tolower= TRUE,
                                         stemming = FALSE,
                                         removeNumbers = TRUE,
                                         removePunctuation = TRUE,
                                         removeStripwhitespace = TRUE,
                                         weighting = weightTf,
                                         bounds=list(global=c(100, Inf))))
dtm
dtm.matrix <- dtm |> as.matrix()


split <-  sample(c("Treinamento",
                   "Validacao",
                   "Teste"),prob=c(0.6,0.2,0.2),
                 size = nrow(dados),
                 replace = TRUE)

# Logistic regression
mq <- glmnet(dtm.matrix[split=="Treinamento",], dados$score_fat[split=="Treinamento"],
             family = "binomial", alpha = 1, lambda = 0)

predito_MMQ <-  predict(mq,
                        newx = dtm.matrix[split=="Teste",],
                        type = "response")


# reliab diagram
p1 <- reliab_plot(predito_MMQ[, 1], dados$score_fat[split == "Teste"],
                  nbins = 50,
                  title = "Logistic Regression")



vc_lasso  <-  cv.glmnet(dtm.matrix[split=="Treinamento",],
                        dados$score_fat[split=="Treinamento"],
                        alpha = 1, family = "binomial")
plot(vc_lasso)

predito_lasso <-  predict(vc_lasso, s = vc_lasso$lambda.min,
                          newx = dtm.matrix[split=="Teste",],
                          type = "response")



p2 <- reliab_plot(predito_lasso[, 1], dados$score_fat[split == "Teste"],
                  nbins = 50,
                  title = "LASSO")


p3 <- data.frame(predito_mq = predito_MMQ[,1]) |>
  ggplot(aes(x = predito_mq)) +
  geom_histogram()


p4 <- data.frame(predito_mq = predito_lasso[,1]) |>
  ggplot(aes(x = predito_mq)) +
  geom_histogram()

# joining graphs
ggpubr::ggarrange(p1, p2, ncol = 2, nrow = 1)
ggsave("reliability_diagrams_amazon.pdf",
       path = "figures",
       width = 8,
       height = 4)


# ECE and MCE
# logistic
ece(predito_MMQ[, 1], 
    dados$score_fat[split == "Teste"], nbins = 50) |>
  round(3)
mce(predito_MMQ[, 1], 
    dados$score_fat[split == "Teste"], nbins = 50) |>
  round(3)

# lasso
ece(predito_lasso[, 1], 
    dados$score_fat[split == "Teste"], nbins = 50) |>
  round(3)
mce(predito_lasso[, 1], 
    dados$score_fat[split == "Teste"], nbins = 50) |>
  round(3)


# Simulated example used in report ----------------------------------------
set.seed(1250)
ntrain <- 5000
nvalid <- 2000
ntest <- 3000
X <- runif(ntrain, min = -1, max = 1)
X_test <- runif(ntest, min = -1, max = 1)
X_valid <- runif(nvalid, min = -1, max = 1)
prob_y <- function(x){x^2}

# train data
sim_train <- data.frame(x = X,
                        p_y = prob_y(X)) |>
  mutate(y = rbinom(length(X), 1, p_y))
sim_train |> glimpse()

# test data
test_data <- data.frame(x = X_test,
                        p_y = prob_y(X_test)) |>
  mutate(y = rbinom(length(X_test), 1, p_y))
test_data |> glimpse()

# validation data
sim_valid <- data.frame(x = X_valid,
                        p_y = prob_y(X_valid)) |>
  mutate(y = rbinom(length(X_valid), 1, p_y))
sim_valid |> glimpse()

# original random forest model
mod_rf <- ranger::ranger(y ~ x, data = sim_train, probability = TRUE)


# real probability versus estimated
x_grid <- seq(-1, 1, 0.01)
data.frame(x = x_grid,
           p_y = prob_y(x_grid)) |>
  ggplot(aes(x = x, y = p_y)) +
  geom_line() +
  theme_bw() +
  labs(y = "P(Y = 1|x)",
       x = "x")
ggsave("prob_sim.pdf",
       path = "figures",
       width = 6,
       height = 4)


plot(x_grid,
     prob_y(x_grid),
     type = "l",
     xlab = "x",
     ylab = "P(y|x)")
lines(x_grid, predict(mod_rf,
                      data.frame(x = x_grid), 
                      type = "response")$predictions[, 2])

pred_rf <- predict(mod_rf,
                   data.frame(x = X_test), 
                   type = "response")$predictions[, 2]

valid_pred <- predict(mod_rf,
                      sim_valid, 
                      type = "response")$predictions[, 2]

# reliability diagram for random forest
reliab_plot(pred_rf, test_data$y)
ggsave("reliab_plot_uncalibrated.pdf",
       path = "figures",
       width = 6,
       height = 4)


# ECE e MCE
ece(pred_rf, test_data$y)
mce(pred_rf, test_data$y)

# histogram binning
hist_bin <- CalibratR:::build_hist_binning(sim_train$y, 
                                           mod_rf$predictions[, 2], 
                                           bins = 50)
pred_bin <- CalibratR:::predict_hist_binning(hist_bin, pred_rf)

p1 <- reliab_plot(pred_bin$predictions, test_data$y, title = "Histogram Binning")
ece(pred_bin$predictions, test_data$y)
mce(pred_bin$predictions, test_data$y)


# isotonic regression
iso_reg <- isoreg(mod_rf$predictions[, 2], sim_train$y) |> as.stepfun()
plot(isoreg(mod_rf$predictions[, 2], sim_train$y))
pred_iso <- iso_reg(pred_rf)

p2 <- reliab_plot(pred_iso, test_data$y, title = "Isotonic Regression")
ece(pred_iso, test_data$y)
mce(pred_iso, test_data$y)


# BBQ
bbq <- CalibratR:::build_BBQ(sim_valid$y, valid_pred)
pred_bbq <- CalibratR:::predict_BBQ(bbq, pred_rf, option = 1)

p3 <- reliab_plot(pred_bbq$predictions, test_data$y, title = "BBQ")
ece(pred_bbq$predictions, test_data$y)
mce(pred_bbq$predictions, test_data$y)

x_grid <- seq(-1, 1, 0.01)
plot(x_grid,
     prob_y(x_grid),
     type = "l",
     xlab = "x",
     ylab = "P(y|x)")
lines(x_grid, predict(mod_rf,
                      data.frame(x = x_grid), 
                      type = "response")$predictions[, 2])
lines(x_grid, CalibratR:::predict_BBQ(bbq, predict(mod_rf,
                                                   data.frame(x = x_grid), 
                                                   type = "response")$predictions[, 2],
                                      option = 1)$predictions, col = "red")

# platt
h <- glm(y ~ g_x, data = sim_valid |>
           mutate(g_x = valid_pred), family = "binomial")

pred_platt <- predict(h,
                      test_data |>
                        mutate(g_x = pred_rf),
                      type = "response")

p4 <- reliab_plot(pred_platt, test_data$y, title = "Platt Scaling")
ece(pred_platt, test_data$y)
mce(pred_platt, test_data$y)


# grouping all measures on a table
ece_s <- map_dbl(list(pred_rf, pred_bin$predictions, pred_iso, pred_bbq$predictions,
                      pred_platt),
                 function(.x){ece(.x, test_data$y, nbins = 50)})

mce_s <- map_dbl(list(pred_rf, pred_bin$predictions, pred_iso, pred_bbq$predictions,
                      pred_platt),
                 function(.x){mce(.x, test_data$y, nbins = 50)})

data_ece <- data.frame(mod = c("uncalibrated",
                               "hist bin",
                               "isoreg",
                               "bbq",
                               "platt")) |>
  mutate(ece = ece_s,
         mce = mce_s)

data_ece |> glimpse()

# joining all reliab diagrams
figure <- ggpubr::ggarrange(p1 + rremove("ylab") + rremove("xlab"), 
                            p2 + rremove("ylab") + rremove("xlab"), 
                            p3 + rremove("ylab") + rremove("xlab"), 
                            p4+ rremove("ylab") + rremove("xlab"),
                            labels = NULL,
                            ncol = 2, nrow = 2,
                            common.legend = TRUE, legend = "bottom",
                            align = "hv", 
                            font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))

annotate_figure(figure, left = textGrob("Calibration Function", rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
                bottom = textGrob("Probability Estimate", gp = gpar(cex = 1.3)))
ggsave("reliab_plots_calibrated.pdf",
       path = "figures",
       width = 10,
       height = 6)
