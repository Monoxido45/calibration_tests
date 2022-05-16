# preliminary results for research project
# same simulation as done before
library(tidyverse)
library(ggpubr)
library(CalibratR)
set.seed(1250)
ntrain <- 5000
nvalid <- 3000
ntest <- 3000
prob_y <- function(x){x^2}

sets <- c("train" = ntrain,
          "valid" = nvalid,
          "test" = ntest) |>
  map(runif, min = -1, max = 1) |>
  map(function(.x){data.frame(x = .x) |>
      mutate(p_y = prob_y(x),
             y = rbinom(.x, 1, p_y))})

train <- sets |> pluck("train")
valid <- sets |> pluck("valid")
test <-  sets |> pluck("test")

# function test if our estimator is A-calibrated

# functions for ECE, MCE and reliab plot
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
  
  plot <- bins_stats |>
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
  
  return(list("bins_stats" = bins_stats,
              "plot" = plot))
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

# brier score
brier_score <- function(prob_estim, ground_truth){
  (prob_estim - ground_truth)^2 |> mean()
}


sim_a_calib <- function(train, valid, test, nbins = 50){
  # training random forest model
  g_x <- ranger::ranger(y ~ x, 
                        data = train,
                        classification = TRUE,
                        probability = TRUE)
  
  # test and validation predictions
  test <- test %>% mutate(g_x = predict(g_x, test,
                                       type = "response")$predictions[, 2],
                          split = case_when(x >= -1 & x < -0.5 ~ "1",
                                            x >= -0.5 & x < 0 ~ "2",
                                            x >= 0 & x < 0.5 ~ "3",
                                            x >= 0.5 & x <= 1 ~ "4"))
    
  
  valid <- valid |> mutate(g_x = predict(g_x, valid,
                                       type = "response")$predictions[, 2],
                           split = case_when(x >= -1 & x < -0.5 ~ "1",
                                             x >= -0.5 & x < 0 ~ "2",
                                             x >= 0 & x < 0.5 ~ "3",
                                             x >= 0.5 & x <= 1 ~ "4"))
  
  # calibrating predictions using the validation set
  # by platt
  h_x <- glm(y ~ g_x, 
             data = valid, family = "binomial")
  
  pred_platt <- predict(h_x, test, type = "response")
  
  # by histogram binning
  hist_bin <- CalibratR:::build_hist_binning(valid$y, 
                                             valid$g_x, 
                                             bins = nbins)
  pred_bin <- CalibratR:::predict_hist_binning(hist_bin, 
                                               test$g_x)$predictions
  
  # plotting all into 4 reliability plot
  # dividing test set A-intervals and then using map
  plot_list <- test |>
    mutate(split = case_when(x >= -1 & x < -0.5 ~ "1",
                             x >= -0.5 & x < 0 ~ "2",
                             x >= 0 & x < 0.5 ~ "3",
                             x >= 0.5 & x <= 1 ~ "4"), 
           pred_bin = pred_bin,
           pred_platt = pred_platt) |>
    group_by(split) |>
    group_split(split) |>
    map_dfr(function(.x){
      p1 <- reliab_plot(.x$g_x, .x$y) |> 
        pluck("bins_stats") |> 
        mutate(mod = "Uncalibrated")
      p2 <- reliab_plot(.x$pred_bin, .x$y) |> 
        pluck("bins_stats") |> 
        mutate(mod = "HB")
      p3 <- reliab_plot(.x$pred_platt, .x$y) |> 
        pluck("bins_stats") |> 
        mutate(mod = "Platt")
      bind_rows(p1, p2 ,p3) |> mutate(partition = unique(.x$split))})
  
  plot_1 <- plot_list |>
    ggplot(aes(x = conf, y = acc)) +
    geom_point(colour = "dodgerblue3",
               alpha = 0.5) +
    geom_line(colour = "dodgerblue3") +
    geom_linerange(aes(x = conf, ymin = conf, 
                       ymax = acc), width = .01)+
    geom_abline(intercept = 0,
                linetype = "dashed",
                colour = "red") +
    theme_bw()+
    scale_x_continuous(breaks = scales::pretty_breaks(8),
                       limits = c(-0.05,1.05))+
    scale_y_continuous(breaks = scales::pretty_breaks(8),
                       limits = c(-0.05,1.05)) +
    facet_grid(row = vars(mod), col = vars(partition)) +
    labs(y = "Calibration Function",
         x = "Probability estimate")
  
  # performing local calibration
  # by platt
  h_x <- glm(y ~ g_x + split, 
             data = valid, 
             family = "binomial")
  
  pred_loc_platt <- predict(h_x, test, type = "response")
  
  # hist bin
  loc_hist <- loc_hist_binning(valid,
                                bins = nbins)
  
  pred_loc_bin <- predict_loc_hist(loc_hist, 
                               test$split,
                               test$g_x)
  
  # plotting reliability plots the same way as before
  # dividing test set A-intervals and then using map
  plot_list <- test |>
    mutate(pred_loc_bin = pred_loc_bin,
           pred_loc_platt = pred_loc_platt) |>
    group_by(split) |>
    group_split(split) |>
    map_dfr(function(.x){
      p1 <- reliab_plot(.x$pred_loc_bin, .x$y) |> 
        pluck("bins_stats") |> 
        mutate(mod = "HB")
      p2 <- reliab_plot(.x$pred_loc_platt, .x$y) |> 
        pluck("bins_stats") |> 
        mutate(mod = "Platt")
      bind_rows(p1, p2) |> mutate(partition = unique(.x$split))})
  
  plot_1 <- plot_list |>
    ggplot(aes(x = conf, y = acc)) +
    geom_point(colour = "dodgerblue3",
               alpha = 0.5) +
    geom_line(colour = "dodgerblue3") +
    geom_linerange(aes(x = conf, ymin = conf, 
                       ymax = acc), width = .01)+
    geom_abline(intercept = 0,
                linetype = "dashed",
                colour = "red") +
    theme_bw()+
    scale_x_continuous(breaks = scales::pretty_breaks(8),
                       limits = c(-0.05,1.05))+
    scale_y_continuous(breaks = scales::pretty_breaks(8),
                       limits = c(-0.05,1.05)) +
    facet_grid(row = vars(mod), col = vars(partition)) +
    labs(y = "Calibration Function",
         x = "Probability estimate")
  
}

loc_hist_binning <- function(valid_data, bins){
  valid_data |> pull(split) |> unique() |> str_sort() |>
  map(function(.x){
      filtered <- valid_data |> 
        filter(split == .x)
      hist_bin <- CalibratR:::build_hist_binning(filtered$y, 
                                     filtered$g_x, 
                                     bins = bins)
        
    })
}


predict_loc_hist <- function(loc_hist_obj, part, preds){
  1:length(preds) |>
    map_dbl(function(.x){
      CalibratR:::predict_hist_binning(loc_hist_obj |> pluck(
        as.numeric(part[.x])),
                                       preds[.x])$predictions
        
    })
}
