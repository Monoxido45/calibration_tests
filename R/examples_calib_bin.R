library(tidyverse)
# biblioteca para alguns metodos de calibracao
library(CalibratR)

# dados simulados
set.seed(1250)
nsim <- 50000
ntest <- 30000
X <- runif(nsim)
X_test <- runif(ntest)
prob_y <- function(x){x^2}

sim_data <- data.frame(x = X,
                       p_y = prob_y(X)) |>
  mutate(y = rbinom(length(X), 1, p_y))
sim_data |> glimpse()

mod_logis <- glm(y ~ x, data = sim_data, family = "binomial")


# formato da probabilidade real versus
x_grid <- seq(0, 1, 0.01)
plot(x_grid,
     prob_y(x_grid),
     type = "l",
     xlab = "x",
     ylab = "P(y|x)")
lines(x_grid, predict(mod_logis,
                      data.frame(x = x_grid), 
                      type = "response"))

# dados do teste
test_data <- data.frame(x = X_test,
                        p_y = prob_y(X_test)) |>
  mutate(y = rbinom(length(X_test), 1, p_y))
test_data |> glimpse()

pred_logis <- predict(mod_logis,
                      test_data,
                      type = "response")


# definindo o grafico de confianca
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

# calculando as medidas para ambos os casos
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

# grafico para o modelo logistico
reliab_plot(pred_logis, test_data$y)


# Exemplo das reviews da amazon -------------------------------------------
set.seed(1)
library(data.table)
library(tm)
library(glmnet)
# por enquanto usando setwd
setwd("~/estatistica_UFSCAR/Doutorado/scripts_test_R")

dados <-  fread("Reviews.csv",
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

# modelo MQ
mq <- glmnet(dtm.matrix[split=="Treinamento",], dados$score_fat[split=="Treinamento"],
             family = "binomial", alpha = 1, lambda = 0)

predito_MMQ <-  predict(mq,
                        newx = dtm.matrix[split=="Teste",],
                        type = "response")


# grafico de confianca
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

# juntando
ggpubr::ggarrange(p1, p2, ncol = 2, nrow = 1)
ggsave("reliability_diagrams_amazon.pdf",
       path = "figures",
       width = 8,
       height = 4)

# calculando medidas
# para logistica
ece(predito_MMQ[, 1], 
    dados$score_fat[split == "Teste"], nbins = 50) |>
  round(3)
mce(predito_MMQ[, 1], 
    dados$score_fat[split == "Teste"], nbins = 50) |>
  round(3)

# para lasso
ece(predito_lasso[, 1], 
    dados$score_fat[split == "Teste"], nbins = 50) |>
  round(3)
mce(predito_lasso[, 1], 
    dados$score_fat[split == "Teste"], nbins = 50) |>
  round(3)







# exemplo simulado --------------------------------------------------------
set.seed(1250)
ntrain <- 5000
nvalid <- 2000
ntest <- 3000
X <- runif(ntrain, min = -1, max = 1)
X_test <- runif(ntest, min = -1, max = 1)
X_valid <- runif(nvalid, min = -1, max = 1)
prob_y <- function(x){x^2}

sim_train <- data.frame(x = X,
                       p_y = prob_y(X)) |>
  mutate(y = rbinom(length(X), 1, p_y))
sim_train |> glimpse()

# dados do teste
test_data <- data.frame(x = X_test,
                        p_y = prob_y(X_test)) |>
  mutate(y = rbinom(length(X_test), 1, p_y))
test_data |> glimpse()

# dados da validacao

sim_valid <- data.frame(x = X_valid,
                        p_y = prob_y(X_valid)) |>
  mutate(y = rbinom(length(X_valid), 1, p_y))
sim_valid |> glimpse()

# modelo original nos dados de treino
mod_rf <- ranger::ranger(y ~ x, data = sim_train, probability = TRUE)


# formato da probabilidade real versus estimado no treino
x_grid <- seq(-1, 1, 0.01)
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

# diagrama de confianca para o modelo nao calibrado
reliab_plot(pred_rf, test_data$y)
ggsave()


# ECE e MCE
ece(pred_rf, test_data$y)
mce(pred_rf, test_data$y)

# RF bem ajustado mas ainda pode ser melhorado
# histogram binning
hist_bin <- CalibratR:::build_hist_binning(sim_train$y, 
                                           mod_rf$predictions[, 2], 
                                           bins = 50)
pred_bin <- CalibratR:::predict_hist_binning(hist_bin, pred_rf)

reliab_plot(pred_bin$predictions, test_data$y)
ece(pred_bin$predictions, test_data$y)
mce(pred_bin$predictions, test_data$y)

# comparando predicoes do calibrado e nao calibrado contra ground truth
x_grid <- seq(-1, 1, 0.01)
plot(x_grid,
     prob_y(x_grid),
     type = "l",
     xlab = "x",
     ylab = "P(y|x)")
lines(x_grid, predict(mod_rf,
                      data.frame(x = x_grid), 
                      type = "response")$predictions[, 2],
      col = "red")
lines(x_grid, CalibratR:::predict_hist_binning(hist_bin, predict(mod_rf,
                                                   data.frame(x = x_grid), 
                    type = "response")$predictions[, 2])$predictions, 
      col = "blue")


# isotonic regression
iso_reg <- isoreg(mod_rf$predictions[, 2], sim_train$y) |> as.stepfun()
pred_iso <- iso_reg(pred_rf)

reliab_plot(pred_iso, test_data$y, nbins = 70)
ece(pred_iso, test_data$y, nbins = 70)
mce(pred_iso, test_data$y, nbins = 70)


# BBQ
bbq <- CalibratR:::build_BBQ(sim_valid$y, valid_pred)
pred_bbq <- CalibratR:::predict_BBQ(bbq, pred_rf, option = 1)

reliab_plot(pred_bbq$predictions, test_data$y)
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

reliab_plot(pred_platt, test_data$y)
ece(pred_platt, test_data$y)
mce(pred_platt, test_data$y)


# todas as medidas em tabela
ece_s <- map_dbl(list(pred_rf, pred_bin$predictions, pred_iso, pred_bbq$predictions,
      pred_platt),
    function(.x){ece(.x, test_data$y, nbins = 70)})

mce_s <- map_dbl(list(pred_rf, pred_bin$predictions, pred_iso, pred_bbq$predictions,
                    pred_platt),
               function(.x){mce(.x, test_data$y, nbins = 70)})

data_ece <- data.frame(mod = c("uncalibrated",
                               "hist bin",
                               "isoreg",
                               "bbq",
                               "platt")) |>
  mutate(ece = ece_s,
         mce = mce_s)

data_ece |> glimpse()
