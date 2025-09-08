# custom anesrake functions 
anesrake_cust <- function (inputter, dataframe, caseid, weightvec = NULL, cap = 5, 
                           verbose = FALSE, maxit = 100, type = "pctlim", pctlim = 5, 
                           nlim = 5, filter = 1, choosemethod = "total", iterate = TRUE, 
                           convcrit = 0.01, force1 = TRUE, center.baseweights = TRUE) 
{
  dataframe <- dataframe[filter == 1, ]
  caseid <- caseid[filter == 1]
  weightvec <- weightvec[filter == 1]
  mat <- as.data.frame(dataframe)
  origtype = type
  fullvars <- 0
  if (is.null(weightvec)) {
    weightvec <- rep(1, length(caseid))
  }
  if (center.baseweights == TRUE) {
    weightvec <- weightvec/mean(weightvec, na.rm = TRUE)
  }
  if (length(weightvec) != length(caseid)) {
    stop("weight vector does not contain the same number of cases as data frame")
  }
  prevec <- weightvec
  not100 <- NULL
  not100 <- names(inputter)[!(sapply(inputter, function(x) sum(x) %in% 
                                       c(1, 100)))]
  # if (!is.null(not100) & force1 == FALSE & length(not100) > 
  #     0) {
  #   warning(paste("Targets for", not100, "do not sum to 100%. Did you make a typo entering the targets?"))
  #   warning(paste("You can force variables to sum to 1 by setting force1 to 'TRUE'"))
  # }
  if (sum(names(inputter) %in% names(dataframe)) != length(names(inputter))) 
    stop(paste("The names of the target variables should match the names of the data frame they are being matched to. The variable(s) -", 
               names(inputter)[!(names(inputter) %in% names(dataframe))], 
               "- were not found in the data frame"))
  if (force1 == TRUE) {
    if (!is.null(not100) & length(not100) > 0) 
      # warning(paste("Targets for", not100, "do not sum to 100%. Adjusting values to total 100%"))
      inputter <- lapply(inputter, function(x) x/sum(x))
  }
  illegalnegs <- sum(unlist(inputter) < 0)
  discrep1 <- anesrakefinder(inputter, dataframe, weightvec, 
                             choosemethod)
  if (type == "nolim") {
    towers <- inputter
  }
  if (type == "pctlim") {
    towers <- selecthighestpcts(discrep1, inputter, pctlim)
  }
  if (type == "nlim") {
    towers <- selectnhighest(discrep1, inputter, nlim)
  }
  if (type == "nmin") {
    towers <- selecthighestpcts(discrep1, inputter, pctlim, 
                                tostop = 0)
    towers2x <- selectnhighest(discrep1, inputter, nlim)
    if (length(towers) > length(towers2x)) {
      type <- "pctlim"
    }
    if (length(towers) < length(towers2x)) {
      towers <- towers2x
    }
  }
  if (type == "nmax") {
    fullvars <- 0
    discrep1 <- anesrakefinder(inputter, dataframe, weightvec, 
                               choosemethod)
    towers <- selecthighestpcts(discrep1, inputter, pctlim)
    towers2x <- selectnhighest(discrep1, inputter, nlim)
    if (length(towers) > length(towers2x)) {
      towers <- towers2x
      fullvars <- 1
    }
  }
  ranweight <- rakelist_cust(towers, mat, caseid, weightvec, cap, 
                             verbose, maxit, convcrit)
  iterations <- ranweight$iterations
  iter1 <- ranweight$iterations
  weightout <- ranweight$weightvec
  if (type == "pctlim" & iterate == TRUE) {
    ww <- 0
    it <- 0
    while (ww < 1) {
      it <- it + 1
      addtotowers <- selecthighestpcts(anesrakefinder(inputter, 
                                                      dataframe, weightout), inputter, pctlim, tostop = 0, 
                                       warn = 0)
      adders <- addtotowers[!(names(addtotowers) %in% names(towers))]
      tow2 <- c(towers, adders)
      towersx <- towers
      towers <- tow2
      if (sum(as.numeric(names(towersx) %in% names(towers))) == 
          length(towers)) {
        ww <- 1
      }
      if (sum(as.numeric(names(towersx) %in% names(towers))) != 
          length(towers)) {
        if (verbose == TRUE) 
          print(paste("Additional variable(s) off after raking, rerunning with variable(s) included"))
        ranweight <- rakelist_cust(towers, mat, caseid, weightvec, 
                                   cap, verbose, maxit, convcrit)
        weightout <- ranweight$weightvec
      }
      if (it > 10) {
        ww <- 1
      }
      iterations <- ranweight$iterations
    }
  }
  if (type == "nmax" & fullvars == 0 & iterate == TRUE) {
    ww <- 0
    it <- 0
    rundiscrep <- discrep1
    discrep2 <- rep(0, length(discrep1))
    while (ww < 1) {
      it <- it + 1
      rundiscrep <- rundiscrep + discrep2
      discrep2 <- anesrakefinder(inputter, dataframe, weightout)
      addtotowers <- selecthighestpcts(discrep2, inputter, 
                                       pctlim, tostop = 0)
      tow2 <- c(towers, addtotowers)
      towersx <- towers
      towers <- unique(tow2)
      names(towers) <- unique(names(tow2))
      if (sum(as.numeric(names(towersx) %in% names(towers))) == 
          length(towers)) {
        ww <- 1
      }
      if (sum(as.numeric(names(towersx) %in% names(towers))) != 
          length(towers)) {
        if (verbose == TRUE) 
          print(paste("Additional variable(s) off after raking, rerunning with variable(s) included"))
        ranweight <- rakelist_cust(towers, mat, caseid, weightvec, 
                                   cap, verbose, maxit, convcrit)
        weightout <- ranweight$weightvec
      }
      if (sum(as.numeric(names(towersx) %in% names(towers))) > 
          nlim) {
        print("variable maximum reached, running on most discrepant overall variables")
        towers <- selectnhighest(discrep1, inputter, 
                                 nlim)
        ranweight <- rakelist_cust(towers, mat, caseid, weightvec, 
                                   cap, verbose, maxit, convcrit)
        weightout <- ranweight$weightvec
        iterations <- 0
        ww <- 1
      }
      if (it >= 10) {
        ww <- 1
      }
      iterations <- ranweight$iterations
    }
  }
  names(weightout) <- caseid
  out <- list(weightvec = weightout, type = type, caseid = caseid, 
              varsused = names(towers), choosemethod = choosemethod, 
              converge = ranweight$converge, nonconvergence = ranweight$nonconvergence, 
              targets = inputter, dataframe = dataframe, iterations = iterations, 
              iterate = iterate, prevec = prevec)
  class(out) <- c("anesrake", "anesrakelist")
  out
}

rakelist_cust <- function (inputter, dataframe, caseid, weightvec = NULL, cap = 999999, 
                           verbose = FALSE, maxit = 1000, convcrit = 0.01) 
{
  mat <- dataframe
  if (is.null(weightvec)) {
    weightvec <- rep(1, length(caseid))
  }
  prevec <- weightvec
  if (sum(is.na(weightvec)) > 0) {
    stop("seed weights cannot have missing values, use filter to eliminate missing values or substitute 1 for missing cases")
  }
  if (length(weightvec) != length(caseid)) {
    stop("weight vector does not contain the same number of cases as data frame")
  }
  if (cap <= 1) {
    stop("cap may not be less than or equal to 1")
  }
  if (cap < 1.5) {
    print("cap is very low, the model may take a long time to run")
  }
  diferr <- 9999999
  diferrold <- 99999999999
  g <- 0
  pctstill <- 1 - convcrit
  pop <- 0
  while (diferr < pctstill * diferrold) {
    g <- g + 1
    wvold <- weightvec
    if (verbose == TRUE) {
      print(paste("Raking...Iteration", g))
    }
    for (i in names(inputter)) {
      weightvec <- rakeonvar(mat[, i], inputter[[i]], weightvec)
    }
    q <- 0
    while (range(weightvec)[2] > cap + 1e-04) {
      q <- q + 1
      if (verbose == TRUE) {
        print(paste("Capping...Iteration ", g, ".", q, 
                    sep = ""))
      }
      weightvec <- sapply(weightvec, function(x) if (x > 
                                                     cap) {
        x <- cap
      }
      else {
        x <- x
      }, simplify = TRUE)
      weightvec <- weightvec/mean(weightvec)
    }
    if (g %in% seq(100, 10000, 50)) {
      print(paste(g, "iterations have occurred, convergence may not be possible...still working"))
    }
    diferrold <- diferr
    diferr <- sum(abs(weightvec - wvold))
    if (verbose == TRUE) {
      print(paste("Current iteration changed total weights by", 
                  diferr))
    }
    if (g > maxit) {
      print(paste("convergence did not occur in", maxit, 
                  "iterations"))
      print("output may not be accurate")
      warning("Raking Algorithm Did Not Converge, Results May Be Highly Inconsistent")
      diferrold <- 0
      pop <- 2
      converge <- paste("No convergence in", maxit, "iterations")
    }
  }
  if (diferr > 0.001) {
    # print("raking achieved only partial convergence, please check the results to ensure that sufficient convergence was achieved.")
    # print(paste("no improvement was apparent after", g, "iterations"))
    # print(paste("current total change in the iteration is:",  diferr, "average change per weight is:", diferr/sum(weightvec)))
    # warning(paste("Raking algorithm achieved only partial convergence, please check the results to ensure that sufficient convergence was achieved.  Average change in weight per case is",  diferr/sum(weightvec)))
    # warning("Results are stable, but do not perfectly match population marginals")
    diferrx <- diferr
    pop <- 1
    converge <- "Results are stable, but do not perfectly match population marginals"
  }
  if (pop == 0) {
    # print(paste("Raking converged in", g, "iterations"))
    diferrx <- diferr
    converge <- "Complete convergence was achieved"
  }
  names(weightvec) <- caseid
  out <- list(weightvec = weightvec, caseid = caseid, iterations = g, 
              nonconvergence = diferr, converge = converge, varsused = names(inputter), 
              targets = inputter, dataframe = dataframe, prevec = prevec)
  class(out) <- "anesrakelist"
  out
}

# raking
require(anesrake)
RAKING <- function(df, reference, GT, sel_col, coalitions, 
                   B = 1, jack = FALSE, sim = FALSE){
  # check right columns selected
  testit::assert(all(sel_col %in% colnames(reference)))
  targets <- NULL
  # change to factors
  # df <- as.data.frame(apply(df, 2, factor))
  df <- df %>% mutate_all(as.factor)
  for (col in sel_col){
    # groupby to get marginal for raking 
    piece <- reference %>% group_by_at(col) %>% summarise('value' = sum(count))
    # make sum to 1
    piece$value <- piece$value / sum(piece$value)
    # add to target vector 
    targets <- c(targets, list(piece$value))
    # add names to values in target 
    names(targets[[length(targets)]]) <- unlist(piece[, 1])
  }
  # assign variable names for pastvote as well 
  names(targets) <- sel_col
  # rake and obtain weights 
  raking_weights <- tryCatch({anesrake_cust(targets, df, 1:nrow(df), cap = 5, pctlim = 0.01, verbose = FALSE)$weightvec}, error = function(e) {
    message("An error occurred in raking: ", e$message)
    raking_weights <- rep(1/nrow(df), nrow(df))
  })
  # obtain weighted mean
  weight_sum <- as.data.frame(questionr::wtd.table(df$nowvote, weights = raking_weights))
  weight_sum$Freq <- weight_sum$Freq / sum(weight_sum$Freq)
  # print results
  raking <- weight_sum
  if ((B > 1) | (jack)){
    if ((B > 1) & (jack)){print('both bootstrap and jacknife indicated, are you sure?')}
    # to append 
    out <- data.frame()
    if (B > 1){to_arrive <- B; print('variance estimation for raking with bootstrap..')}
    if (jack){to_arrive <- nrow(df); print('variance estimation for raking with jacknife..')}
    # for (i in 1:nrow(df)){
    for (i in 1:to_arrive){
      cat("\rProcessing", i, " our of ", nrow(df))  # \r moves the cursor to the beginning of the line
      flush.console()
      if (B > 1){
        jack_idx <- sample(1:nrow(df), nrow(df), replace = TRUE)
        jack_df <- df[jack_idx, ]
      }
      if (jack){jack_df <- df[c(1:(i-1), (i+1):nrow(df)), ]}
      # rake and obtain weights 
      suppressWarnings(
        raking_weights <- anesrake_cust(targets, jack_df, 1:nrow(jack_df), cap = 5, pctlim = 0.01, verbose = FALSE)$weightvec
      )
      # obtain weighted mean 
      weight_sum <- as.data.frame(questionr::wtd.table(jack_df$nowvote, weights = raking_weights))
      weight_sum$Freq <- weight_sum$Freq / sum(weight_sum$Freq)
      # generate to append 
      tmp <- weight_sum$Freq
      names(tmp) <- weight_sum$Var1
      out <- rbind(out, t(as.data.frame(tmp, stringsAsFactors = FALSE)))
    }
  }
  # obtain standard error
  raking$se <- rep(NA, nrow(raking))
  if (B > 1){raking$se <- apply(out, 2, sd)}
  if (jack){raking$se <- ((nrow(df) - 1) / nrow(df)) * sqrt(apply(sweep(out, 2, raking$Freq, "-")^2, 2, sum))}
  raking$A <- raking$Freq
  raking$Freq <- NULL
  rownames(raking) <- raking$Var1 
  raking$Var1 <- NULL
  # if ()
  # raking[coalitions, c('A', 'se')]
  result_processor(raking[coalitions, c('A', 'se')], coalitions, GT, df, B, jack) 
}

# IPW 
make_reference_sample <- function(poptable, sample_size = 20000){
  library(tidyr)
  # tmp_census <- poptable$census %>% group_by_at(sel_col) %>% summarise('count' = sum(count))
  poptable$n <- round(poptable$count * sample_size, 0)
  # Use uncount() to repeat rows based on 'value' column
  reference_sample <- (poptable %>% uncount(n))
  reference_sample
}

PSIPW <- function(df, reference_sample, GT, sel_col, coalitions, B = 1, jack = FALSE){
  # concatenate non-probability sample with reference sample
  concatenated_df <- bind_rows(df[sel_col], reference_sample[sel_col])
  # add weights 
  concatenated_df$weights <- c(rep(1, nrow(df[sel_col])), reference_sample$count)
  # add indicator variable 
  concatenated_df$S <- c(rep(1, nrow(df[sel_col])), rep(0, nrow(reference_sample)))
  # adjust weights for size of non-probability dataset
  concatenated_df$weights <- concatenated_df$weights * ((sum(1 / reference_sample$count) - nrow(df)) / sum(1 / reference_sample$count))
  # estimate model  # CORRECT HERE MODEL FORMULA 
  ipw_formula <- as.formula(paste0('S ~ ', paste(sel_col, collapse = ' + ')))
  mod_ipw <- glm(ipw_formula, data = concatenated_df, weights = weights,
                 family = binomial(link = 'logit')) # add weights 
  # obtain weights (inverse of propensity)
  ipw_w <- 1 / predict(mod_ipw, newdata = df, type = 'response')
  # obtain weighted mean 
  weight_sum <- as.data.frame(questionr::wtd.table(df$nowvote, weights = ipw_w[1:nrow(df)]))
  weight_sum$Freq <- weight_sum$Freq / sum(weight_sum$Freq)
  ipw <- weight_sum
  if ((B > 1) | (jack)){
    if ((B > 1) & (jack)){print('both bootstrap and jacknife indicated, are you sure?')}
    out <- data.frame()
    if (B > 1){to_arrive <- B; print('variance estimation for PSIPW with bootstrap..')}
    if (jack){to_arrive <- nrow(df); print('variance estimation for PSIPW with jacknife..')}
    for (i in 1:to_arrive){
      cat("\rProcessing", i, " our of ", nrow(df))  # \r moves the cursor to the beginning of the line
      flush.console()
      if (B > 1){
        jack_idx <- sample(1:nrow(df), nrow(df), replace = TRUE)
        jack_df <- df[jack_idx, ]
      }
      if (jack){jack_df <- df[c(1:(i-1), (i+1):nrow(df)), ]}
      m <- nrow(jack_df[sel_col])
      # concatenate non-probability sample with reference sample
      concatenated_jack <- bind_rows(jack_df[sel_col], reference_sample[sel_col])
      # add weights 
      concatenated_jack$weights <- c(rep(1, m), reference_sample$count)
      # add indicator variable 
      concatenated_jack$S <- c(rep(1, m), rep(0, nrow(reference_sample)))
      # adjust weights for size of non-probability dataset
      concatenated_jack$weights <- concatenated_jack$weights * ((sum(1 / reference_sample$count) - m) / sum(1 / reference_sample$count))
      # estimate model 
      suppressWarnings(
        mod_ipw <- glm(ipw_formula, data = concatenated_jack, weights = weights, family = binomial(link = 'logit')) # add weights 
      )
      # obtain weights (inverse of propensity)
      ipw_w <- 1 / predict(mod_ipw, newdata = jack_df, type = 'response')
      # obtain weighted mean 
      weight_sum <- as.data.frame(questionr::wtd.table(jack_df$nowvote, weights = ipw_w[1:m]))
      weight_sum$Freq <- weight_sum$Freq / sum(weight_sum$Freq)
      # generate to append 
      tmp <- weight_sum$Freq
      names(tmp) <- weight_sum$Var1
      out <- rbind(out, t(as.data.frame(tmp, stringsAsFactors = FALSE)))
    }
  }
  # obtain standard error
  ipw$se <- rep(NA, nrow(ipw))
  if (B > 1){ipw$se <- apply(out, 2, sd)}
  if (jack){ipw$se <- ((nrow(df) - 1) / nrow(df)) * sqrt(apply(sweep(out, 2, ipw$Freq, "-")^2, 2, sum))}
  ipw$A <- ipw$Freq
  rownames(ipw) <- ipw$Var1 
  result_processor(ipw[coalitions, c('A', 'se')], coalitions, GT, df, B = B, jack = jack)
}

# post-strat 
POSTSTRAT <- function(df, reference, GT, sel_col, coalitions, B = 1){
  reference <- reference %>% group_by_at(sel_col) %>% summarise(count = sum(count))
  poststrat_formula <- paste0('as.factor(nowvote) ~ ', paste(sel_col, collapse = ' + '))
  # multinomial model 
  mod_mult <- nnet::multinom(poststrat_formula, data = df, trace = FALSE)
  # post-strat
  post_strat <- predict(mod_mult, newdata = reference, type = "probs")
  poststrat <- apply(post_strat * reference$count / sum(reference$count), 2, sum)
  poststrat <- as.data.frame(poststrat)
  poststrat$A <- poststrat$poststrat
  out <- data.frame()
  # bootstrap variance estimate
  if (B > 1){
    for (b in 1:B){
      # cat("\rProcessing", b, " our of ", B)  # \r moves the cursor
      # flush.console()
      boot_index <- sample(nrow(df), nrow(df), replace = TRUE)
      boot_df <- df[boot_index, ]
      suppressMessages(mod_mult <- nnet::multinom(poststrat_formula, data = boot_df, trace = FALSE))
      post_strat <- predict(mod_mult, newdata = reference, type = "probs")
      tmp <- apply(post_strat * reference$count / sum(reference$count), 2, sum)
      out <- rbind(out, t(as.data.frame(tmp, stringsAsFactors = FALSE)))
    }
    # add standard error 
    poststrat$se <- apply(out, 2, sd)
  } else {poststrat$se <- rep(NA, nrow(poststrat))}
  result_processor(poststrat[coalitions, ], coalitions, GT, df, B = B, jack = FALSE)
}

# MRP
MRP <- function(df, reference, GT, sel_col, coalitions, B = 1, method = 'lme4', 
                chains = 4, iter = 2000, cores = 1, refresh = 0){
  if ((B > 1) & (method == 'rstanarm')){print('bootstrap indicated but bayesian estimate')}
  mrp <- NULL
  for (party in coalitions){
    # print(party)
    MPR_formula <- as.formula(paste0('nowvote == party ~ ', paste(paste0('(1|', sel_col, ')'), collapse = ' + ')))
    if (method == 'rstanarm'){
      require(rstanarm)
      # estimate bayesian model
      suppressWarnings(
        MRP <- rstanarm::stan_glmer(MPR_formula, data = df, family = binomial(link = 'logit'),
                                    chains = chains, iter = iter, refresh = refresh, cores = cores)
      )
      # posterior predict 
      epred_mat <- posterior_epred(MRP, newdata = reference, draws = 50)
      # post-stratification step
      mrp_estimates_vector <- epred_mat %*% reference$count / sum(reference$count)
      mrp[[party]] <- c(A = mean(mrp_estimates_vector), se = sd(mrp_estimates_vector))
    }
    if (method == 'lme4'){
      require(lme4)
      # print(party)
      MPR_formula <- as.formula(paste0('nowvote == party ~ ', paste(paste0('(1|', sel_col, ')'), collapse = ' + ')))
      MRP <- lme4::glmer(MPR_formula, data = df, family = binomial(link = 'logit'))
      pred_helper <- function(MRP){
        sum(predict(MRP, newdata = reference, type = 'response') * reference$count)
      }
      boot_est <- lme4::bootMer(MRP, pred_helper, nsim = B)
      mrp[[party]] <- c(A = boot_est$t0, se = sd(boot_est$t))
    }
  }
  result_processor(as.data.frame(t(as.data.frame(mrp)))[coalitions, ],
                   coalitions, GT, df, jack = TRUE)
}

# DRP
DRP <- function(df, reference, GT, sel_col, coalitions){
  require(survey)
  require(nonprobsvy)
  reference <- reference %>% group_by_at(sel_col) %>% summarise(count = sum(count))
  # pass data to binomial dummy format
  df_dummy <- fastDummies::dummy_cols(df[sel_col], remove_first_dummy = TRUE,
                                      remove_selected_columns = TRUE)
  # same for census
  census_dummy <- fastDummies::dummy_cols(reference, 
                                          select_columns = sel_col,
                                          remove_first_dummy = TRUE, 
                                          remove_selected_columns = TRUE)
  freq_dummy <- apply(census_dummy, 2, function(x) sum(x * census_dummy$count))
  # N of voters in the reference population
  N <- 45210950
  pop_totals <- as.integer(freq_dummy * N)
  names(pop_totals) <- c('Y', paste0('X', seq(1, length(names(freq_dummy))-1)))
  # set pop totals variable for synthjoin 
  pop_totals['`(Intercept)`'] <- N
  colnames(df_dummy) <- paste0('X', seq(1, length(colnames(df_dummy))))
  out <- NULL
  for (party in coalitions){
    df_dummy$Y <- as.integer(df$nowvote == party)
    mod_col <- paste0('X', seq(1, length(colnames(df_dummy))-1))
    # Estimate DRP
    suppressWarnings(
      dr_logit_poptotals <- nonprob(
        data = df_dummy,
        selection = as.formula(paste0('~ ', paste(mod_col, collapse = ' + '))),
        outcome = as.formula(paste0('Y ~ ', paste(mod_col, collapse = ' + '))),
        pop_totals = pop_totals[c(length(pop_totals), 2:(length(pop_totals)-1))],
        method_selection = "logit",
        # svydesign = sample_prob, # for reference sample case
        # control_selection = controlSel(est_method_sel = "gee", h = 1)
      )
    )
    # obtain point estimate and variance
    out[[party]] <- cbind(dr_logit_poptotals$output,dr_logit_poptotals$confidence_interval)
  }
  drp <- bind_rows(out)
  rownames(drp) <- coalitions
  drp$A <- drp$mean
  drp$se <- drp$SE
  # jack = true just to use se, not actual Jacknife
  result_processor(drp[coalitions, c('A', 'se')], coalitions, GT, df, jack = TRUE)
}

# ---- Neural Net ---
torch_NN_mn <- function(ref,                  # reference census 
                        GT,                   # ground truth
                        db,                   # dataframe to perform analysis on
                        cols,                 # columns to select
                        coalitions,           # political coalitions
                        training_size = 0.9,  # training size
                        max_patience = 3,     # training epochs of patience
                        tolerance = 0.0005,   # tolerance in validation loss change
                        n_epochs = 20,        # max number of epochs
                        batch_size = 12,      # batch size
                        learning_rate = 1e-2, # learning rate
                        n_neurons = 4,        # n of neurons in each layer 
                        loss_type = 'new',    # type of loss function
                        early_stopping = TRUE # activate early stopping
){
  sq_loss <- function(y, X){
    (X - y)^2
  }

  # select data without target
  net_db <- as.data.frame(db[, cols])
  # categorical to dummy for neural net
  x_ <- as.data.frame(fastDummies::dummy_cols(net_db) %>% select_if(~ all(. %in% 0:1)))
  # select train 
  train_indices <- sample(c(TRUE, FALSE), nrow(x_), 
                          replace=TRUE, prob=c(training_size, 1 - training_size))
  # get train and test
  x_train <- as.matrix(x_[train_indices, ])
  x_test <- as.matrix(x_[!train_indices, ])
  # pass to tensor 
  x_torch <- torch_tensor(x_train, dtype=torch_float32())
  x_torch_v <- torch_tensor(x_test, dtype=torch_float32())
  # define input size 
  n_len <- ncol(x_train)
  # define output size
  n_out <- length(coalitions)
  # re-factor in right order
  to_dummy <- factor(db$nowvote, levels = coalitions)
  # change to dummy 
  y_ <- fastDummies::dummy_cols(to_dummy) %>% select_if(~ all(. %in% 0:1))
  # pass to tensor 
  y_train <- as.matrix(y_[train_indices, ])
  y_test <- as.matrix(y_[!train_indices, ])
  y_torch <- torch_tensor(y_train, dtype=torch_float32())
  y_torch_v <- torch_tensor(y_test, dtype=torch_float32())
  # define nets
  h_net <- nn_sequential(
    nn_linear(n_len, n_neurons), nn_relu(), 
    nn_linear(n_neurons, n_neurons), nn_relu(), 
    nn_linear(n_neurons, n_out), nn_softmax(-1) 
  )
  # init weights using Xavier method 
  init_xavier(h_net)
  init_xavier(a_net)
  # set parameters and optimizer
  params <- c(h_net$parameters)
  # optimizer
  optimizer_ru <- optim_adam(params, lr = learning_rate)
  # define control variables 
  losses_log <- NULL; outer_loop <- FALSE; patience <- 0; 
  # start training
  for (epoch in seq(n_epochs)){
    for (i in seq(0, nrow(x_train), batch_size)){
      Xbatch <- x_torch[i:(i + batch_size), ]
      h <- h_net(Xbatch)
      a <- a_net(Xbatch)
      ybatch <- y_torch[i : (i + batch_size)]
      epoch_loss <- sq_loss(ybatch, h, a, directions, gammas) 
      if ((epoch == 1) & (i == 0)){
        # print(sprintf('Initial loss for party %s is %g', sel_coal, ru_loss))
        # get initial validation loss 
        h_v <- h_net(x_torch_v)
        val_loss <- sq_loss(ybatch, h)
        losses_log <- c(losses_log, as.numeric(val_loss))
        last_loss <- val_loss
      }
      optimizer_ru$zero_grad()
      epoch_loss$backward()
      optimizer_ru$step()
    }
    # get validation loss
    h_v <- h_net(x_torch_v)
    val_loss <- sq_loss(ybatch, h)
    if (early_stopping){
      val_loss <- sq_loss(y_torch_v, h_v)
    }
    # print(sprintf("Finished epoch %g, validation loss is %g", epoch, val_loss))
    # check for over-training
    if ((as.numeric(abs(val_loss - last_loss)) < tolerance) | 
        (as.numeric(val_loss) > as.numeric(last_loss))){
      patience <- patience + 1
      # print(sprintf('patience is now %g', patience))
      if (patience >= max_patience){
        # print(sprintf('interrupting..'))
        outer_loop <- TRUE
        break 
      }
    } else {patience <- 0}
    losses_log <- c(losses_log, as.numeric(val_loss))
    last_loss <- val_loss
    if (outer_loop) break
  }
  # adjust reference 
  ref <- ref %>% group_by_at(cols) %>% summarise(value = sum(count))
  # obtain dummy for census
  to_pred <- as.matrix(fastDummies::dummy_cols(ref[, cols[cols %in% colnames(ref)]]) %>%
                         select_if(~ all(. %in% 0:1)))
  # check columns are in order 
  to_pred <- to_pred[, colnames(x_train)]
  testit::assert(all(colnames(x_train) == colnames(to_pred)))
  # output predicted post-stratified quantity 
  to_pred <- h_net(torch_tensor(to_pred, dtype=torch_float32()))
  # get prediction
  out <- to_pred * torch_tensor(ref$value, 
                                dtype=torch_float32())$reshape(c(-1, 1))
  out <- apply(out, 2, sum)
  names(out) <- coalitions
  dru <- as.data.frame(out)
  dru$se <- NA
  dru$A <- dru$out
  result_processor(dru[coalitions, ], coalitions, GT, db)
}
