# pequeno teste para o doutorado
library(tidyverse)
# teste doq a gnt quer fazer comparado com oq tem
# Função que roda simulação -----------------------------------------------
# roda simulacao
# adicionando reliability diagram
roda_sim <- function(X, X_test, prob_y, pred_real, nbins = 50, seed = 1250){
  set.seed(seed)
  set.seed(125)
  # dados da simulacao
  sim_data <- data.frame(x = X,
                         p_y = prob_y(X)) |>
    mutate(y = rbinom(length(X), 1, p_y))
  
  # dados do teste
  test_data <- data.frame(x = X_test,
                          p_y = prob_y(X_test)) |>
    mutate(y = rbinom(length(X_test), 1, p_y))
  
  # treino e validacao
  split <- sim_data %>% rsample::initial_split(strata = "y", prop = 0.7)
  
  sim_train <- split %>% rsample::training()
  sim_valid <- split %>% rsample::testing()
  
  
  # modelo g_x
  if(exists("g_x")){rm(g_x)}
  g_x <- ranger::ranger(y ~ x, 
                        data = sim_train,
                        classification = TRUE,
                        probability = TRUE)
  
  
  # predizendo com g(x) e x como covariavel
  sim_valid <- sim_valid %>%
    mutate(g_x = predict(g_x, sim_valid,
                         type = "response")$predictions[, 2])
  
  h_x <- glm(y ~ x + g_x, 
             data = sim_valid, family = "binomial")
  
  # ajustando g(x) usando um knn
  x_grid <- seq(0, 1, 0.001)
  
  
  pred_h <- function(x_test){
    pred_g_x <- predict(g_x, data.frame(x = x_test),
                        type = "response")$predictions[, 2]
    
    predict(h_x, newdata = data.frame(x = x_test,
                                      g_x = pred_g_x),
            type = "response")
    
  }
  
  pred_g <- function(x_test){
    predict(g_x, data.frame(x = x_test),
            type = "response")$predictions[, 2]
  }
  
  # h so com g_x
  
  pred_h_g_x <- function(x_test){
    pred_g_x <- predict(g_x, data.frame(x = x_test),
                        type = "response")$predictions[, 2]
    h_x <- glm(y ~ g_x, 
               data = sim_valid, family = "binomial")
    
    predict(h_x, newdata = data.frame(g_x = pred_g_x),
            type = "response")
    
  }
  
  pred_f <- function(x_test){
    h_x <- glm(y ~ x, 
               data = sim_valid, family = "binomial")
    
    predict(h_x, newdata = data.frame(x = x_test),
            type = "response")
  }
  
  
  p1 <- data.frame(x = x_grid) %>%
    ggplot(aes(x = x))+
    stat_function(fun = pred_g, aes(colour = "g"))+
    stat_function(fun = pred_h, aes(colour = "h(x, g(x))"))+
    stat_function(fun = pred_h_g_x, aes(colour = "h(g(x))"))+
    stat_function(fun = pred_f, aes(colour = "f"))+
    stat_function(fun = pred_real, aes(colour = "real"))+
    theme_bw() +
    labs(y = "Probabilidade estimada",
         x = "x",
         colour = "Ajustes") +
    scale_colour_brewer(palette = "Set1")
  
  p2 <- h_x |> broom::tidy() %>%
    filter(term != "(Intercept)") %>%
    ggplot(aes(x = term, y = estimate)) +
      geom_bar(fill = "dodgerblue3", colour = "dodgerblue3", 
               stat = "identity", alpha = 0.65) +
      theme_bw() +
      labs(y = "Valor",
           x = "Coeficiente")
  
  # diagramas de confianca
  x_grid <- seq(0, 1, length.out = nbins)
  b_m <- x_grid[2] - x_grid[1]
  # criando as estatisticas dos bins para o h_g_x
  bins_stats <- x_grid[-1] |>
    map_dfr(~ test_data |>
          mutate(pred_bin_h = pred_h(x)) |>
          filter(pred_bin_h <= .x & pred_bin_h > (.x - b_m))|>
          summarise(acc = mean((y == 1)),
                    conf = mean(pred_bin_h)) |>
          mutate(gap = abs(acc - conf)) |>
          mutate(lim_sup =.x))
  p3 <- bins_stats |>
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
  
  # criando as estatisticas dos bins para o h_x
  bins_stats <- x_grid[-1] |>
    map_dfr(~ test_data |>
              mutate(pred_bin_h_g = pred_h_g_x(x)) |>
              filter(pred_bin_h_g <= .x & pred_bin_h_g > (.x - b_m))|>
              summarise(acc = mean((y == 1)),
                        conf = mean(pred_bin_h_g)) |>
              mutate(gap = abs(acc - conf)) |>
              mutate(lim_sup =.x))

  
  p4 <- bins_stats |>
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
         title = "h(g(x))")+
    scale_x_continuous(breaks = scales::pretty_breaks(8),
                       limits = c(-0.05,1.05))+
    scale_y_continuous(breaks = scales::pretty_breaks(8),
                       limits = c(-0.05,1.05))
    
  ggpubr::ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
}


# experimento com X^2 -----------------------------------------------------
# gerando os dados
set.seed(1250)
nsim <- 5000
ntest <- 3000
X <- runif(nsim)
X_test <- runif(ntest)
prob_y <- function(x){x^2}
pred_real <- function(x_test){
  x_test^2
}

roda_sim(X, X_test, prob_y, pred_real)



# experimento com X^3 -----------------------------------------------------
set.seed(1250)
nsim <- 5000
ntest <- 3000
X <- runif(nsim)
X_test <- runif(ntest)
prob_y <- function(x){x^3}
pred_real <- function(x_test){
  x_test^3
}

roda_sim(X, X_test, prob_y, pred_real)

# experimento com X^5 -----------------------------------------------------
set.seed(1250)
nsim <- 5000
ntest <- 3000
X <- runif(nsim)
X_test <- runif(ntest)
prob_y <- function(x){x^5}
pred_real <- function(x_test){
  x_test^5
}

roda_sim(X, X_test, prob_y, pred_real)



