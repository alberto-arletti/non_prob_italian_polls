# import required packages
library(dplyr)
library(stats)
library(torch)
library(survey)

# prepare shared post-stratification inputs
prepare_poststrat_inputs <- function(data, census,
                                     outcome_var = "nowvote",
                                     count_var = "count") {
  # check that the census count column exists
  if (!(count_var %in% names(census))) stop("Missing in census: ", count_var)
  
  # infer post-stratification variables from census
  strata_vars <- setdiff(names(census), count_var)
  
  # stop if no post-stratification variables are available
  if (length(strata_vars) == 0) stop("No post-stratification variables found in census.")
  
  # collect required data columns
  data_cols <- c(strata_vars, outcome_var)
  
  # check missing columns in survey data
  missing_data_cols <- setdiff(data_cols, names(data))
  
  # stop if survey data lacks census variables or outcome
  if (length(missing_data_cols) > 0) stop("Missing in data: ", paste(missing_data_cols, collapse = ", "))
  
  # create helper for full-cell identifiers
  make_postcell <- function(df, vars) {
    do.call(interaction, c(df[vars], list(drop = TRUE, sep = "___"))) |> as.character()
  }
  
  # remove incomplete survey rows
  data0 <- data |>
    # keep complete cases on needed survey variables
    filter(if_all(all_of(data_cols), ~ !is.na(.x)))
  
  # remove incomplete census rows
  census0 <- census |>
    # keep complete cases on needed census variables
    filter(if_all(all_of(c(strata_vars, count_var)), ~ !is.na(.x)))
  
  # harmonise survey variables
  data0 <- data0 |>
    # avoid factor-level mismatch
    mutate(across(all_of(strata_vars), as.character))
  
  # harmonise census variables
  census0 <- census0 |>
    # avoid factor-level mismatch
    mutate(across(all_of(strata_vars), as.character))
  
  # prepare population cells
  pop_df <- census0 |>
    # keep post-stratification variables and count
    select(all_of(c(strata_vars, count_var))) |>
    # aggregate duplicate census cells
    group_by(across(all_of(strata_vars))) |>
    # sum population counts
    summarise(pop_count = sum(.data[[count_var]], na.rm = TRUE), .groups = "drop")
  
  # create population cell id
  pop_df$postcell <- make_postcell(pop_df, strata_vars)
  
  # create survey cell id
  data0$postcell <- make_postcell(data0, strata_vars)
  
  # keep only positive population cells
  pop_positive <- pop_df |>
    # remove zero-population cells
    filter(pop_count > 0)
  
  # remove survey rows in absent or zero-population cells
  data_use <- data0 |>
    # retain only valid population cells
    semi_join(pop_positive |> select(postcell), by = "postcell")
  
  # stop if no rows are left
  if (nrow(data_use) == 0) stop("No survey rows remain after matching to positive census cells.")
  
  # keep only population cells represented in the survey
  pop_use <- pop_positive |>
    # retain sample-covered cells
    filter(postcell %in% data_use$postcell)
  
  # rescale population totals to retained sample size
  pop_use <- pop_use |>
    # create frequency totals for survey functions
    mutate(Freq = pop_count / sum(pop_count) * nrow(data_use))
  
  # define common postcell levels
  common_levels <- sort(unique(pop_use$postcell))
  
  # align survey postcell factor
  data_use$postcell <- factor(data_use$postcell, levels = common_levels)
  
  # align population postcell factor
  pop_use$postcell <- factor(pop_use$postcell, levels = common_levels)
  
  # make outcome categorical
  data_use[[outcome_var]] <- factor(data_use[[outcome_var]])
  
  # return prepared objects
  list(data = data_use, population = pop_use)
}

# format survey estimates
format_svy_estimates <- function(est, outcome_var = "nowvote",
                                 conf_level = 0.95) {
  # compute confidence intervals
  ci <- confint(est, level = conf_level)
  
  # return clean estimate table
  data.frame(
    # outcome category
    category = sub(paste0("^", outcome_var), "", names(coef(est))),
    # adjusted estimate
    estimate = as.numeric(coef(est)),
    # standard error
    se = as.numeric(SE(est)),
    # confidence interval lower bound
    ci_low = ci[, 1],
    # confidence interval upper bound
    ci_high = ci[, 2],
    # clean row names
    row.names = NULL
  )
}

# post-stratification using survey::postStratify()
poststrat <- function(data, census,
                      outcome_var = "nowvote",
                      count_var = "count",
                      conf_level = 0.95) {
  # prepare common inputs
  prep <- prepare_poststrat_inputs(data, census, outcome_var, count_var)
  
  # create base survey design
  des0 <- svydesign(
    # no clustering
    ids = ~1,
    # equal base weights
    weights = ~1,
    # filtered survey data
    data = prep$data
  )
  
  # run full post-stratification
  des_ps <- postStratify(
    # base design
    design = des0,
    # full-cell post-stratum
    strata = ~postcell,
    # known population totals by cell
    population = prep$population |> select(postcell, Freq)
  )
  
  # estimate outcome shares
  est <- svymean(as.formula(paste0("~", outcome_var)), des_ps)
  
  # return estimates only
  format_svy_estimates(est, outcome_var, conf_level)
}

# calibration using survey::calibrate()
calib <- function(data, census,
                  outcome_var = "nowvote",
                  count_var = "count",
                  conf_level = 0.95,
                  calfun = "linear",
                  bounds = c(0.01, 100)) {
  # prepare common inputs
  prep <- prepare_poststrat_inputs(data, census, outcome_var, count_var)
  
  # create base survey design
  des0 <- svydesign(
    # no clustering
    ids = ~1,
    # equal base weights
    weights = ~1,
    # filtered survey data
    data = prep$data
  )
  
  # create calibration model matrix
  mm <- model.matrix(~postcell - 1, data = prep$data)
  
  # create population totals for postcell dummies
  pop_totals <- prep$population$Freq
  
  # name totals to match model matrix columns
  names(pop_totals) <- paste0("postcell", prep$population$postcell)
  
  # reorder population totals to match model matrix
  pop_totals <- pop_totals[colnames(mm)]
  
  # stop if names failed to match
  if (any(is.na(pop_totals))) stop("Population totals do not match calibration matrix columns.")
  
  # check finite bounds for logit calibration
  if (calfun == "logit" && (length(bounds) != 2 || any(!is.finite(bounds)))) {
    stop("Logit calibration requires finite bounds, e.g. bounds = c(0.01, 100).")
  }
  
  # check valid bound order
  if (calfun == "logit" && bounds[1] >= bounds[2]) {
    stop("For logit calibration, bounds[1] must be smaller than bounds[2].")
  }
  
  # run calibration with finite bounds if logit is selected
  if (calfun == "logit") {
    des_cal <- calibrate(
      # base design
      design = des0,
      # one dummy per post-stratification cell
      formula = ~postcell - 1,
      # known cell totals
      population = pop_totals,
      # bounded calibration distance
      calfun = "logit",
      # finite weight bounds
      bounds = bounds
    )
  } else {
    des_cal <- calibrate(
      # base design
      design = des0,
      # one dummy per post-stratification cell
      formula = ~postcell - 1,
      # known cell totals
      population = pop_totals,
      # calibration distance function
      calfun = calfun
    )
  }
  
  # estimate outcome shares
  est <- svymean(as.formula(paste0("~", outcome_var)), des_cal)
  
  # return estimates only
  format_svy_estimates(est, outcome_var, conf_level)
}


library(nonprobsvy)

# generate approximate reference sample from census
generate_reference <- function(census, nrow = 20000,
                               count_var = "count",
                               seed = NULL) {
  # set seed if requested
  if (!is.null(seed)) set.seed(seed)
  
  # check count column
  if (!(count_var %in% names(census))) stop("Missing in census: ", count_var)
  
  # infer X variables from census
  x_vars <- setdiff(names(census), count_var)
  
  # stop if no X variables are available
  if (length(x_vars) == 0) stop("No reference variables found in census.")
  
  # clean census table
  census0 <- census |>
    # keep complete rows
    filter(if_all(all_of(c(x_vars, count_var)), ~ !is.na(.x))) |>
    # keep positive counts
    filter(.data[[count_var]] > 0) |>
    # harmonise variables
    mutate(across(all_of(x_vars), as.character))
  
  # create sampling probabilities
  prob <- census0[[count_var]] / sum(census0[[count_var]])
  
  # sample census cells approximately
  idx <- sample(seq_len(base::nrow(census0)), size = nrow, replace = TRUE, prob = prob)
  
  # expand sampled cells into reference survey
  reference <- census0[idx, x_vars, drop = FALSE]
  
  # add equal reference weights
  reference$ref_weight <- 1
  
  # clean row names
  rownames(reference) <- NULL
  
  # return reference sample
  reference
}

# extract one numeric vector from nonprobsvy output
extract_nonprob_vector <- function(fit, what, target_names) {
  # extract requested object
  tmp <- as.data.frame(extract(fit, what))
  
  # return named columns if present
  if (all(target_names %in% names(tmp))) return(as.numeric(tmp[1, target_names]))
  
  # return named rows if present
  if (all(target_names %in% rownames(tmp))) {
    # find numeric columns
    num_cols <- names(tmp)[sapply(tmp, is.numeric)]
    
    # return first numeric column
    return(as.numeric(tmp[target_names, num_cols[1]]))
  }
  
  # return first numeric column if order matches
  if (nrow(tmp) == length(target_names)) {
    # find numeric columns
    num_cols <- names(tmp)[sapply(tmp, is.numeric)]
    
    # return first numeric column
    return(as.numeric(tmp[[num_cols[1]]]))
  }
  
  # stop if format is unexpected
  stop("Could not extract ", what, " from nonprobsvy output.")
}

# pooled-sample IPW with bootstrap SEs and optional weight capping
ps_ipw <- function(df, reference,
                   outcome_var = "nowvote",   # name of the Y / party column in df
                   sel_col = NULL,            # common predictors; NULL = auto-infer
                   count_var = "count",       # population count column in reference
                   trim = 5,                  # cap weights at this value T; NULL = no cap
                   B = 200,                   # number of bootstrap replicates
                   conf_level = 0.95,         # for the normal-approx CI
                   return_weights = FALSE) {  # attach 1st-rep uncapped weights as attribute
  # auto-infer predictors: columns shared by both, minus outcome and count
  if (is.null(sel_col)) {
    sel_col <- setdiff(intersect(names(df), names(reference)), c(outcome_var, count_var))
  }
  # reference gets unit counts if none supplied
  ref_w <- if (count_var %in% names(reference)) reference[[count_var]] else rep(1, nrow(reference))
  # keep complete rows in the non-prob sample
  df <- df[stats::complete.cases(df[c(sel_col, outcome_var)]), , drop = FALSE]
  # keep complete rows in the reference
  ok <- stats::complete.cases(reference[sel_col]) & !is.na(ref_w)
  reference <- reference[ok, , drop = FALSE]; ref_w <- ref_w[ok]
  # harmonise each predictor to a factor with levels shared across both samples
  for (v in sel_col) {
    lev <- sort(unique(c(as.character(df[[v]]), as.character(reference[[v]]))))
    df[[v]] <- factor(as.character(df[[v]]), levels = lev)
    reference[[v]] <- factor(as.character(reference[[v]]), levels = lev)
  }
  # fix the party set once so every replicate has the same columns
  parties <- sort(unique(as.character(df[[outcome_var]])))
  # outcome as character for fast comparison
  y_all <- as.character(df[[outcome_var]])
  # reference block of the pooled model frame (constant across replicates)
  ref_block <- reference[sel_col]
  # selection-model formula: S ~ v1 + v2 + ...
  f <- stats::reformulate(sel_col, response = "S")
  # one replicate: bootstrap df, fit logit, weight, return party proportions
  one_rep <- function(idx, keep_w = FALSE) {
    # bootstrap rows of the non-prob sample
    dfb <- df[idx, sel_col, drop = FALSE]; yb <- y_all[idx]
    # pooled predictors: non-prob on top, reference below
    pooled <- rbind(dfb, ref_block)
    # indicator: 1 = non-prob, 0 = reference
    pooled$S <- c(rep(1, length(idx)), rep(0, nrow(ref_block)))
    # case weights as a column (1 for non-prob, population counts for reference)
    pooled$.cw <- c(rep(1, length(idx)), ref_w)
    # fit the propensity model (suppress the non-integer-weight binomial warning)
    mod <- suppressWarnings(stats::glm(f, data = pooled, weights = .cw,
                                       family = binomial(link = "logit")))
    # predicted P(S = 1 | X) on the bootstrap non-prob rows
    p <- stats::predict(mod, newdata = dfb, type = "response")
    # IPW weight = 1 / propensity (as in your example)
    w <- 1 / p
    # normalise to mean 1 so the cap T means "multiples of the average weight"
    w <- w / mean(w)
    # stash the uncapped (mean-1) weights if requested
    if (keep_w) wenv$first_w <- as.numeric(w)
    # cap weights at T if trimming is on
    if (!is.null(trim)) w <- pmin(w, trim)
    # weighted proportion per fixed party (0 if a party is absent this draw)
    vapply(parties, function(pt) sum(w * (yb == pt)) / sum(w), numeric(1))
  }
  # environment to hold the uncapped weights of the first replicate
  wenv <- new.env()
  # matrix of bootstrap estimates: parties x B
  reps <- vapply(seq_len(B), function(b) {
    # resample non-prob row indices with replacement
    idx <- sample.int(nrow(df), replace = TRUE)
    # compute proportions; capture weights only on the first replicate
    one_rep(idx, keep_w = (b == 1 && return_weights))
  }, numeric(length(parties)))
  # bootstrap mean estimate per party (sums to 1)
  estimate <- rowMeans(reps)
  # bootstrap SE per party
  se <- apply(reps, 1, stats::sd)
  # normal critical value
  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  # tidy output table
  out <- data.frame(category = parties,
                    estimate = as.numeric(estimate),
                    se = as.numeric(se),
                    ci_low = as.numeric(estimate) - z * se,
                    ci_high = as.numeric(estimate) + z * se,
                    row.names = NULL, stringsAsFactors = FALSE)
  # optionally attach the first replicate's uncapped weight distribution
  if (return_weights) {out <- list(out, wenv$first_w)}
  # return the tidy table
  out
}

# doubly robust regression and post-stratification with nonprobsvy
nonprob_drp <- function(data, census,
                        outcome_var = "nowvote",
                        count_var = "count",
                        N = 45210950,
                        method_selection = "logit",
                        conf_level = 0.95,
                        return_fit = FALSE) {
  # reuse helper_rev.R preparation
  prep <- prepare_poststrat_inputs(data, census, outcome_var, count_var)
  
  # use all census covariates in both models
  x_vars <- setdiff(names(census), count_var)
  
  # keep only variables used by nonprobsvy
  data0 <- prep$data |>
    dplyr::select(dplyr::all_of(c(outcome_var, x_vars)))
  
  # keep retained population cells
  pop0 <- prep$population
  
  # define the selection model
  selection_formula <- stats::reformulate(x_vars)
  
  # define the binary outcome model
  outcome_formula <- stats::reformulate(x_vars, response = "Y")
  
  # build population model matrix
  X_pop <- stats::model.matrix(selection_formula, data = pop0)
  
  # rescale census counts/proportions to voter population size
  pop_n <- pop0$pop_count / sum(pop0$pop_count) * N
  
  # compute population totals for nonprobsvy
  pop_totals <- colSums(X_pop * pop_n)
  
  # request analytical standard errors
  inf_control <- control_inf(
    var_method = "analytic",
    alpha = 1 - conf_level
  )
  
  # get outcome categories
  outcome_levels <- levels(data0[[outcome_var]])
  
  # confidence interval multiplier
  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  
  # fit one binary DRP model per category
  fit_one <- function(lvl) {
    # copy data for this category
    data_loop <- data0
    
    # create binary outcome
    data_loop$Y <- as.integer(data_loop[[outcome_var]] == lvl)
    
    # estimate DRP
    fit <- suppressWarnings(nonprob(
      data = data_loop,
      selection = selection_formula,
      outcome = outcome_formula,
      pop_totals = pop_totals,
      method_selection = method_selection,
      control_inference = inf_control,
      se = TRUE
    ))
    
    # extract adjusted estimate
    est <- extract_nonprob_vector(fit, "mean", "Y")
    
    # extract analytical standard error
    se <- extract_nonprob_vector(fit, "se", "Y")
    
    # return row and fitted object
    list(
      row = data.frame(
        category = lvl,
        estimate = est,
        se = se,
        ci_low = est - z * se,
        ci_high = est + z * se,
        row.names = NULL
      ),
      fit = fit
    )
  }
  
  # run all binary models
  res <- lapply(outcome_levels, fit_one)
  
  # bind rows into helper_rev.R output shape
  out <- dplyr::bind_rows(lapply(res, `[[`, "row"))
  
  # return fits only when requested
  if (return_fit) {
    fits <- lapply(res, `[[`, "fit")
    names(fits) <- outcome_levels
    return(list(estimates = out, fits = fits))
  }
  
  # return estimates only
  out
}

# unweighted raw vote-share estimate with multinomial uncertainty
uniweighted <- function(data, outcome_var = "nowvote",
                        conf_level = 0.95,
                        ci_method = c("wilson", "wald"),
                        categories = NULL) {
  # match confidence interval method
  ci_method <- match.arg(ci_method)
  
  # use now_vote if requested default is absent
  if (!outcome_var %in% names(data) && outcome_var == "nowvote" && "now_vote" %in% names(data)) outcome_var <- "now_vote"
  
  # stop if outcome column is missing
  if (!outcome_var %in% names(data)) stop("Outcome variable not found: ", outcome_var)
  
  # keep original outcome vector
  y_raw <- data[[outcome_var]]
  
  # drop missing outcomes
  y <- y_raw[!is.na(y_raw)]
  
  # stop if no valid outcomes remain
  if (length(y) == 0) stop("No non-missing outcomes found in: ", outcome_var)
  
  # preserve supplied categories or factor levels
  if (is.null(categories) && is.factor(y_raw)) categories <- levels(y_raw)
  
  # otherwise infer observed categories
  if (is.null(categories)) categories <- sort(unique(as.character(y)))
  
  # count all requested categories, including zero-count ones
  tab <- table(factor(as.character(y), levels = categories))
  
  # store effective sample size for the outcome
  n <- sum(tab)
  
  # convert counts to shares
  p_hat <- as.numeric(tab) / n
  
  # analytic multinomial standard errors
  se <- sqrt(p_hat * (1 - p_hat) / n)
  
  # compute normal quantile
  z <- qnorm(1 - (1 - conf_level) / 2)
  
  # create Wald intervals
  if (ci_method == "wald") {
    ci_low <- pmax(0, p_hat - z * se)
    ci_high <- pmin(1, p_hat + z * se)
  }
  
  # create Wilson score intervals
  if (ci_method == "wilson") {
    denom <- 1 + z^2 / n
    centre <- (p_hat + z^2 / (2 * n)) / denom
    half_width <- z * sqrt((p_hat * (1 - p_hat) / n) + z^2 / (4 * n^2)) / denom
    ci_low <- pmax(0, centre - half_width)
    ci_high <- pmin(1, centre + half_width)
  }
  
  # build multinomial covariance matrix
  vcov_mat <- (diag(p_hat) - tcrossprod(p_hat)) / n
  
  # name covariance rows
  rownames(vcov_mat) <- names(tab)
  
  # name covariance columns
  colnames(vcov_mat) <- names(tab)
  
  # return clean estimate table
  out <- data.frame(
    # outcome category
    category = names(tab),
    # number of observations
    # n = as.numeric(tab),
    # unweighted estimate
    estimate = p_hat,
    # multinomial standard error
    se = se,
    # lower confidence interval
    ci_low = ci_low,
    # upper confidence interval
    ci_high = ci_high,
    # clean row names
    row.names = NULL
  )
  
  # attach full covariance matrix for later aggregate metrics
  attr(out, "vcov") <- vcov_mat
  
  # return estimate table
  out
}

# multilevel regression and post-stratification using rstanarm
mrp <- function(data, census,
                outcome_var = "nowvote",
                count_var = "count",
                draws = 50,
                chains = 4,
                iter = 2000,
                cores = 1,
                refresh = 0,
                conf_level = 0.95,
                return_fit = FALSE) {
  # check that rstanarm is available
  if (!requireNamespace("rstanarm", quietly = TRUE)) stop("Install rstanarm first.")
  # reuse helper_rev.R cleaning and support checks
  prep <- prepare_poststrat_inputs(data, census, outcome_var, count_var)
  # infer adjustment variables from the population table
  x_vars <- setdiff(names(census), count_var)
  # keep only model variables in the sample
  data0 <- prep$data |> dplyr::select(dplyr::all_of(c(outcome_var, x_vars)))
  # prepare full population table for post-stratification
  pop0 <- census |>
    # keep positive and complete population cells
    dplyr::filter(dplyr::if_all(dplyr::all_of(c(x_vars, count_var)), ~ !is.na(.x)), .data[[count_var]] > 0) |>
    # harmonise predictors before aggregation
    dplyr::mutate(dplyr::across(dplyr::all_of(x_vars), as.character)) |>
    # aggregate duplicate population cells
    dplyr::group_by(dplyr::across(dplyr::all_of(x_vars))) |>
    # sum population counts
    dplyr::summarise(pop_count = sum(.data[[count_var]], na.rm = TRUE), .groups = "drop")
  # align sample and population factor levels
  for (v in x_vars) {
    # collect shared modelling levels
    lev <- sort(unique(c(data0[[v]], pop0[[v]])))
    # apply levels to sample and population
    data0[[v]] <- factor(data0[[v]], levels = lev); pop0[[v]] <- factor(pop0[[v]], levels = lev)
  }
  # make the outcome categorical
  data0[[outcome_var]] <- factor(data0[[outcome_var]])
  # store population weights for post-stratification
  pop_w <- pop0$pop_count / sum(pop0$pop_count)
  # build random-intercept formula
  mrp_formula <- stats::as.formula(paste0("Y ~ ", paste(paste0("(1|", x_vars, ")"), collapse = " + ")))
  # fit one binary model per outcome category
  fit_one <- function(lvl) {
    # create binary outcome for this category
    data_loop <- dplyr::mutate(data0, Y = as.integer(.data[[outcome_var]] == lvl))
    # fit Bayesian multilevel logistic regression
    fit <- suppressWarnings(rstanarm::stan_glmer(mrp_formula, 
                                                 data = data_loop, 
                                                 family = stats::binomial("logit"), 
                                                 chains = chains, 
                                                 iter = iter, 
                                                 cores = cores, 
                                                 refresh = refresh))
    # predict posterior probabilities for population cells
    epred <- rstanarm::posterior_epred(fit, newdata = pop0, draws = draws, allow_new_levels = TRUE)
    # post-stratify posterior draws
    theta <- as.vector(epred %*% pop_w)
    # return row and fitted model
    list(row = data.frame(category = lvl, 
                          estimate = mean(theta), 
                          se = stats::sd(theta), 
                          ci_low = stats::quantile(theta, (1 - conf_level) / 2, names = FALSE), 
                          ci_high = stats::quantile(theta, 1 - (1 - conf_level) / 2, names = FALSE), 
                          row.names = NULL), 
         fit = fit)
  }
  # run all category models
  res <- lapply(levels(data0[[outcome_var]]), fit_one)
  # bind rows into the helper_rev.R estimate shape
  out <- dplyr::bind_rows(lapply(res, `[[`, "row"))
  # optionally return fitted stan_glmer models
  if (return_fit) return(list(estimates = out, 
                              fits = stats::setNames(lapply(res, `[[`, "fit"), 
                                                     levels(data0[[outcome_var]])), 
                              data = data0, population = pop0))
  # return estimates only
  out
}

# model-assisted post-stratification with multinomial logit
multinomial <- function(data, census,
                        outcome_var = "nowvote",
                        count_var = "count",
                        x_vars = NULL,
                        B = 1,
                        conf_level = 0.95) {
  # check nnet availability
  if (!requireNamespace("nnet", quietly = TRUE)) stop("Install nnet first.")
  
  # infer adjustment variables from census
  if (is.null(x_vars)) x_vars <- setdiff(names(census), count_var)
  
  # check required columns
  if (length(setdiff(c(x_vars, outcome_var), names(data))) > 0) stop("Missing columns in data.")
  if (length(setdiff(c(x_vars, count_var), names(census))) > 0) stop("Missing columns in census.")
  
  # clean survey data
  data0 <- data |>
    # keep complete model rows
    filter(if_all(all_of(c(x_vars, outcome_var)), ~ !is.na(.x))) |>
    # harmonise predictors
    mutate(across(all_of(x_vars), as.character))
  
  # clean and collapse population table
  pop0 <- census |>
    # keep positive complete population cells
    filter(if_all(all_of(c(x_vars, count_var)), ~ !is.na(.x)), .data[[count_var]] > 0) |>
    # harmonise predictors
    mutate(across(all_of(x_vars), as.character)) |>
    # collapse duplicate cells
    group_by(across(all_of(x_vars))) |>
    # sum counts
    summarise(pop_count = sum(.data[[count_var]]), .groups = "drop")
  
  # align factor support
  for (v in x_vars) {
    # keep levels available in both data sources
    lev <- sort(intersect(unique(data0[[v]]), unique(pop0[[v]])))
    # restrict to shared levels
    data0 <- data0[data0[[v]] %in% lev, , drop = FALSE]; pop0 <- pop0[pop0[[v]] %in% lev, , drop = FALSE]
    # set common factor levels
    data0[[v]] <- factor(data0[[v]], levels = lev); pop0[[v]] <- factor(pop0[[v]], levels = lev)
  }
  
  # make outcome categorical
  data0[[outcome_var]] <- factor(data0[[outcome_var]])
  
  # store outcome levels
  y_levels <- levels(data0[[outcome_var]])
  
  # create multinomial main-effects formula
  f <- reformulate(x_vars, response = outcome_var)
  
  # helper: fit model and post-stratify
  est_once <- function(dat) {
    # preserve all outcome levels in bootstrap samples
    dat[[outcome_var]] <- factor(dat[[outcome_var]], levels = y_levels)
    # fit multinomial linear model
    mod <- suppressMessages(nnet::multinom(f, data = dat, trace = FALSE))
    # predict population-cell probabilities
    pred <- predict(mod, newdata = pop0, type = "probs")
    # handle binary outcome returned as vector
    if (is.null(dim(pred))) pred <- cbind(1 - pred, pred)
    # enforce outcome names
    colnames(pred) <- y_levels
    # post-stratify over population cells
    colSums(pred[, y_levels, drop = FALSE] * pop0$pop_count / sum(pop0$pop_count))
  }
  
  # point estimate
  est <- est_once(data0)
  
  # bootstrap standard errors
  se <- if (B > 1) {
    # bootstrap repeated estimates
    boot <- replicate(B, est_once(data0[sample(nrow(data0), nrow(data0), replace = TRUE), , drop = FALSE]))
    # category-level standard errors
    apply(boot, 1, sd, na.rm = TRUE)
  } else {
    # no bootstrap requested
    rep(NA_real_, length(est))
  }
  
  # confidence multiplier
  z <- qnorm(1 - (1 - conf_level) / 2)
  
  # return helper_rev.R estimate shape
  data.frame(
    category = names(est),
    estimate = as.numeric(est),
    se = as.numeric(se),
    ci_low = as.numeric(est - z * se),
    ci_high = as.numeric(est + z * se),
    row.names = NULL
  )
}

# one-hot encode categorical predictors for torch models
make_nn_dummy_matrices <- function(data, population, x_vars) {
  # align predictor levels across sample and population
  for (v in x_vars) {
    # collect all observed levels
    lev <- sort(unique(c(as.character(data[[v]]), as.character(population[[v]]))))
    # factorise the sample column
    data[[v]] <- factor(as.character(data[[v]]), levels = lev)
    # factorise the population column
    population[[v]] <- factor(as.character(population[[v]]), levels = lev)
  }
  
  # helper to encode one dataframe
  encode_one <- function(df) {
    # build one dummy block for each predictor
    blocks <- lapply(x_vars, function(v) {
      # create dummy columns for the current predictor
      mm <- stats::model.matrix(~ x - 1, data = data.frame(x = df[[v]]))
      # use readable and stable column names
      colnames(mm) <- paste0(v, "__", sub("^x", "", colnames(mm)))
      # return the dummy block
      mm
    })
    # combine all dummy blocks
    do.call(cbind, blocks)
  }
  
  # encode sample rows
  x_data <- encode_one(data)
  # encode population cells
  x_pop <- encode_one(population)
  # check that columns are aligned
  if (!identical(colnames(x_data), colnames(x_pop))) stop("Sample and population dummy matrices are not aligned.")
  # return encoded objects
  list(x_data = x_data, x_pop = x_pop, data = data, population = population)
}

# initialise linear layers with Xavier weights
init_torch_xavier <- function(net) {
  # loop over sequential modules
  for (i in seq_len(length(net))) {
    # initialise only linear layers
    if (inherits(net[[i]], "nn_linear")) {
      # initialise layer weights
      torch::nn_init_xavier_uniform_(net[[i]]$weight)
      # initialise layer bias
      torch::nn_init_constant_(net[[i]]$bias, 0)
    }
  }
  # return the network invisibly
  invisible(net)
}

# # categorical cross-entropy for one-hot outcomes
# nn_onehot_loss <- function(pred, y, eps = 1e-7) {
#   # avoid log of zero
#   pred <- torch::torch_clamp(pred, min = eps, max = 1 - eps)
#   # compute mean negative log likelihood
#   -torch::torch_mean(torch::torch_sum(y * torch::torch_log(pred), dim = 2))
# }

# categorical cross-entropy for one-hot outcomes
nn_onehot_loss <- function(pred, y, eps = 1e-7) {
  # avoid log of zero
  pred <- torch::torch_clamp(pred, min = eps, max = 1 - eps)
  # mean negative log-likelihood; sum over the class dim (last dim, robust to 1-row batches)
  -torch::torch_mean(torch::torch_sum(y * torch::torch_log(pred), dim = -1))
}

# fit a small torch network and post-stratify predictions
fit_nn_poststrat_once <- function(data, population, x_data, x_pop,
                                  outcome_var, y_levels,
                                  training_size = 0.8,
                                  max_patience = 3,
                                  tolerance = 0.0005,
                                  n_epochs = 50,
                                  batch_size = 32,
                                  learning_rate = 1e-3,
                                  n_neurons = 8,
                                  n_hidden = 2,
                                  dropout = 0,
                                  weight_decay = 0,
                                  early_stopping = TRUE,
                                  verbose = FALSE) {
  # create the one-hot outcome matrix
  y_mat <- stats::model.matrix(~ y - 1, data = data.frame(y = factor(data[[outcome_var]], levels = y_levels)))
  # make sure outcome columns follow y_levels
  colnames(y_mat) <- y_levels
  
  # define train rows
  train_idx <- sample(c(TRUE, FALSE), nrow(x_data), replace = TRUE, prob = c(training_size, 1 - training_size))
  # avoid an empty training split
  if (sum(train_idx) == 0) train_idx[sample(seq_len(nrow(x_data)), 1)] <- TRUE
  # avoid an empty validation split
  if (sum(!train_idx) == 0) train_idx[sample(which(train_idx), 1)] <- FALSE
  
  # create training tensors
  x_train <- torch::torch_tensor(x_data[train_idx, , drop = FALSE], dtype = torch::torch_float32())
  # create validation tensors
  x_valid <- torch::torch_tensor(x_data[!train_idx, , drop = FALSE], dtype = torch::torch_float32())
  # create training outcome tensor
  y_train <- torch::torch_tensor(y_mat[train_idx, , drop = FALSE], dtype = torch::torch_float32())
  # create validation outcome tensor
  y_valid <- torch::torch_tensor(y_mat[!train_idx, , drop = FALSE], dtype = torch::torch_float32())
  
  # start network layer list
  layers <- list(torch::nn_linear(ncol(x_data), n_neurons), torch::nn_relu())
  # add additional hidden layers
  if (n_hidden > 1) {
    for (j in seq_len(n_hidden - 1)) {
      # add optional dropout before the next layer
      if (dropout > 0) layers <- c(layers, list(torch::nn_dropout(p = dropout)))
      # add hidden linear and activation layers
      layers <- c(layers, list(torch::nn_linear(n_neurons, n_neurons), torch::nn_relu()))
    }
  }
  # add output layer and softmax
  layers <- c(layers, list(torch::nn_linear(n_neurons, length(y_levels)), torch::nn_softmax(dim = -1)))
  # create sequential network
  net <- do.call(torch::nn_sequential, layers)
  # initialise weights
  init_torch_xavier(net)
  
  # create optimiser
  opt <- torch::optim_adam(net$parameters, lr = learning_rate, weight_decay = weight_decay)
  # initialise early-stopping trackers
  last_loss <- Inf; patience <- 0; losses <- numeric(0)
  
  # train across epochs
  for (epoch in seq_len(n_epochs)) {
    # set network to training mode
    net$train()
    # loop over mini-batches
    for (i in seq(1, nrow(x_data[train_idx, , drop = FALSE]), by = batch_size)) {
      # define current batch rows
      batch_rows <- i:min(i + batch_size - 1, nrow(x_data[train_idx, , drop = FALSE]))
      # predict current batch
      pred <- net(x_train[batch_rows, ])
      # compute current loss
      loss <- nn_onehot_loss(pred, y_train[batch_rows, ])
      # reset gradients
      opt$zero_grad()
      # back-propagate loss
      loss$backward()
      # update parameters
      opt$step()
    }
    
    # set network to evaluation mode
    net$eval()
    # compute validation predictions
    val_pred <- net(x_valid)
    # compute validation loss
    val_loss <- nn_onehot_loss(val_pred, y_valid)
    # store validation loss
    losses <- c(losses, as.numeric(val_loss))
    
    # print progress if requested
    if (verbose) message("epoch ", epoch, ": validation loss = ", round(as.numeric(val_loss), 6))
    # stop early if validation loss stalls or worsens
    if (early_stopping) {
      # update patience when loss does not improve enough
      if (abs(as.numeric(last_loss) - as.numeric(val_loss)) < tolerance || as.numeric(val_loss) > as.numeric(last_loss)) patience <- patience + 1 else patience <- 0
      # interrupt training once patience is exhausted
      if (patience >= max_patience) break
    }
    # update previous loss
    last_loss <- val_loss
  }
  
  # predict population-cell probabilities
  pred_pop <- net(torch::torch_tensor(x_pop, dtype = torch::torch_float32()))
  # convert predictions to an R matrix
  pred_pop <- as.matrix(as.array(pred_pop))
  # guard binary or single-row shape
  if (is.null(dim(pred_pop))) pred_pop <- matrix(pred_pop, ncol = length(y_levels))
  # enforce outcome names
  colnames(pred_pop) <- y_levels
  
  # compute population weights
  pop_w <- population$pop_count / sum(population$pop_count)
  # post-stratify predicted probabilities
  estimate <- colSums(pred_pop[, y_levels, drop = FALSE] * pop_w)
  # return fit and estimates
  list(estimate = estimate, fit = net, losses = losses)
}

# infer neural-network predictors from shared census and sample columns
infer_nn_x_vars <- function(data, census,
                            outcome_var = "nowvote",
                            count_var = "count") {
  # keep variables available in both files
  x_vars <- intersect(names(census), names(data))
  
  # remove non-predictor columns
  x_vars <- setdiff(x_vars, c(count_var, outcome_var, "postcell", "pop_count", "Freq"))
  
  # stop if nothing can be used
  if (length(x_vars) == 0) stop("No shared census/data predictors found for neural-network post-stratification.")
  
  # return inferred variables
  x_vars
}

nn_poststrat <- function(data, census,
                         outcome_var = "nowvote",
                         count_var = "count",
                         x_vars = NULL,
                         B = 1,
                         conf_level = 0.95,
                         training_size = 0.8,
                         max_patience = 3,
                         tolerance = 0.0005,
                         n_epochs = 50,
                         batch_size = 32,
                         learning_rate = 1e-3,
                         n_neurons = 8,
                         n_hidden = 2,
                         dropout = 0,
                         weight_decay = 0,
                         early_stopping = TRUE,
                         seed = NULL,
                         verbose = FALSE,
                         return_fit = FALSE) {
  # check torch availability
  if (!requireNamespace("torch", quietly = TRUE)) stop("Install torch first.")
  
  # set R and torch seeds if requested
  if (!is.null(seed)) { set.seed(seed); torch_manual_seed(seed) }
  
  # infer adjustment variables from shared census/data columns
  if (is.null(x_vars)) x_vars <- infer_nn_x_vars(data, census, outcome_var, count_var)
  
  # stop if the outcome is missing
  if (!(outcome_var %in% names(data))) stop("Missing in data: ", outcome_var)
  
  # stop if the count is missing
  if (!(count_var %in% names(census))) stop("Missing in census: ", count_var)
  
  # check required survey columns
  missing_data_cols <- setdiff(c(x_vars, outcome_var), names(data))
  
  # stop if survey columns are missing
  if (length(missing_data_cols) > 0) stop("Missing columns in data: ", paste(missing_data_cols, collapse = ", "))
  
  # check required census columns
  missing_census_cols <- setdiff(c(x_vars, count_var), names(census))
  
  # stop if census columns are missing
  if (length(missing_census_cols) > 0) stop("Missing columns in census: ", paste(missing_census_cols, collapse = ", "))
  
  # clean survey data
  data0 <- data |>
    # keep complete model rows
    filter(if_all(all_of(c(x_vars, outcome_var)), ~ !is.na(.x))) |>
    # harmonise predictors
    mutate(across(all_of(x_vars), as.character))
  
  # clean and collapse population table
  pop0 <- census |>
    # keep positive complete population cells
    filter(if_all(all_of(c(x_vars, count_var)), ~ !is.na(.x)), .data[[count_var]] > 0) |>
    # harmonise predictors
    mutate(across(all_of(x_vars), as.character)) |>
    # collapse duplicate cells
    group_by(across(all_of(x_vars))) |>
    # sum population counts
    summarise(pop_count = sum(.data[[count_var]], na.rm = TRUE), .groups = "drop")
  
  # stop if no complete survey rows remain
  if (nrow(data0) == 0) stop("No complete survey rows available for neural-network post-stratification.")
  
  # stop if no population cells remain
  if (nrow(pop0) == 0) stop("No positive population cells available for neural-network post-stratification.")
  
  # keep sample rows inside population support
  for (v in x_vars) {
    # retain only sample levels present in the population table
    data0 <- data0[data0[[v]] %in% unique(pop0[[v]]), , drop = FALSE]
  }
  
  # stop if support matching removed everything
  if (nrow(data0) == 0) stop("No survey rows remain after matching population support.")
  
  # reset row names before bootstrap operations
  rownames(data0) <- NULL
  
  # encode predictors for torch
  enc <- make_nn_dummy_matrices(data0, pop0, x_vars)
  
  # update cleaned sample
  data0 <- enc$data
  
  # update cleaned population
  pop0 <- enc$population
  
  # make outcome categorical
  data0[[outcome_var]] <- factor(data0[[outcome_var]])
  
  # collect outcome levels
  y_levels <- levels(data0[[outcome_var]])
  
  # helper for one complete estimate
  est_once <- function(dat, x_mat) {
    # preserve outcome support in bootstrap samples
    dat[[outcome_var]] <- factor(dat[[outcome_var]], levels = y_levels)
    
    # fit the torch network and post-stratify
    fit_nn_poststrat_once(dat, pop0, x_mat, enc$x_pop,
                          outcome_var = outcome_var,
                          y_levels = y_levels,
                          training_size = training_size,
                          max_patience = max_patience,
                          tolerance = tolerance,
                          n_epochs = n_epochs,
                          batch_size = batch_size,
                          learning_rate = learning_rate,
                          n_neurons = n_neurons,
                          n_hidden = n_hidden,
                          dropout = dropout,
                          weight_decay = weight_decay,
                          early_stopping = early_stopping,
                          verbose = verbose)
  }
  
  # compute point estimate, guarded so a failed full-data fit degrades instead of aborting (highest level)
  point <- tryCatch(est_once(data0, enc$x_data),
                    # on failure fall back to NA estimates over the known outcome levels
                    error = function(e) list(estimate = setNames(rep(NA_real_, length(y_levels)), y_levels), fit = NULL, losses = NULL))
  
  # extract point estimate vector
  est <- point$estimate
  
  # bootstrap standard errors if requested
  se <- if (B > 1) {
    # run bootstrap refits
    boot <- replicate(B, {
      # wrap the entire repetition so one failed refit is skipped, not fatal (highest level)
      tryCatch({
        # sample row positions with replacement
        idx <- sample(seq_len(nrow(data0)), nrow(data0), replace = TRUE)
        
        # create bootstrap data by position
        boot_dat <- data0[idx, , drop = FALSE]
        
        # discard duplicated row names from resampling
        rownames(boot_dat) <- NULL
        
        # create matching bootstrap design matrix by position
        boot_x <- enc$x_data[idx, , drop = FALSE]
        
        # return bootstrap estimate
        est_once(boot_dat, boot_x)$estimate
      },
      # on any error drop this repetition by returning NAs of the right length
      error = function(e) rep(NA_real_, length(est)))
    })
    
    # compute category-level standard errors
    apply(boot, 1, sd, na.rm = TRUE)
  } else {
    # no bootstrap requested
    rep(NA_real_, length(est))
  }
  
  # compute confidence multiplier
  z <- qnorm(1 - (1 - conf_level) / 2)
  
  # build helper_rev.R estimate table
  out <- data.frame(category = names(est),
                    estimate = as.numeric(est),
                    se = as.numeric(se),
                    ci_low = as.numeric(est - z * se),
                    ci_high = as.numeric(est + z * se),
                    row.names = NULL)
  
  # return fit and diagnostics if requested
  if (return_fit) return(list(estimates = out, fit = point$fit, losses = point$losses, data = data0, population = pop0, x_vars = x_vars))
  
  # return estimates only
  out
}

# add to the imports at the top of helper_rev.R
library(xgboost)

# fit one S~X gradient-boosted propensity model and return non-prob propensities
fit_ipw_xgb_once <- function(x_np, x_ref, ref_w,
                             eta = 0.1,
                             max_depth = 4,
                             subsample = 0.8,
                             colsample_bytree = 0.8,
                             min_child_weight = 1,
                             nrounds = 500,
                             early_stopping_rounds = 20,
                             training_size = 0.8,
                             seed = NULL) {
  # pool non-prob (label 1) and reference (label 0) design matrices
  X <- rbind(x_np, x_ref)
  # build the selection indicator
  Z <- c(rep(1, nrow(x_np)), rep(0, nrow(x_ref)))
  # weight reference rows by their survey weight, non-prob rows by 1
  w <- c(rep(1, nrow(x_np)), ref_w)
  # draw a fresh booster seed from the (already seeded) R stream
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1)
  # split pooled rows into train/validation for early stopping
  train_idx <- sample(c(TRUE, FALSE), nrow(X), replace = TRUE, prob = c(training_size, 1 - training_size))
  # guard against an empty training split
  if (!any(train_idx)) train_idx[sample(seq_len(nrow(X)), 1)] <- TRUE
  # guard against an empty validation split
  if (all(train_idx)) train_idx[sample(which(train_idx), 1)] <- FALSE
  # training matrix for xgboost
  dtrain <- xgb.DMatrix(data = X[train_idx, , drop = FALSE], label = Z[train_idx], weight = w[train_idx])
  # validation matrix used to monitor logloss
  dvalid <- xgb.DMatrix(data = X[!train_idx, , drop = FALSE], label = Z[!train_idx], weight = w[!train_idx])
  # boosting hyper-parameters for binary propensity
  params <- list(objective = "binary:logistic", eval_metric = "logloss",
                 eta = eta, max_depth = max_depth, subsample = subsample,
                 colsample_bytree = colsample_bytree, min_child_weight = min_child_weight, seed = seed, nthread = 1)
  # shared training arguments
  args <- list(params = params, data = dtrain, nrounds = nrounds,
               early_stopping_rounds = early_stopping_rounds, verbose = 0)
  # train with early stopping, supporting both new (evals) and old (watchlist) xgboost
  bst <- tryCatch(do.call(xgb.train, c(args, list(evals = list(valid = dvalid)))),
                  error = function(e) do.call(xgb.train, c(args, list(watchlist = list(valid = dvalid)))))
  # use the early-stopped best iteration when predicting, with version-safe fallbacks
  bi <- bst$best_iteration
  # predict propensity p = P(non-prob | x) for the non-prob rows
  p <- if (!is.null(bi)) {
    tryCatch(predict(bst, x_np, iterationrange = c(1, bi + 1)),
             error = function(e) tryCatch(predict(bst, x_np, ntreelimit = bst$best_ntreelimit),
                                          error = function(e2) predict(bst, x_np)))
  } else predict(bst, x_np)
  # return propensities and fitted model
  list(p = as.numeric(p), fit = bst)
}

# inverse-probability weighting with a gradient-boosted (XGBoost) selection model
nonprob_ipw_xgb <- function(data, reference_survey,
                            outcome_var = "nowvote",
                            selection_vars = NULL,
                            ref_weight_var = "ref_weight",
                            conf_level = 0.95,
                            var_method = "bootstrap",
                            num_boot = 500,
                            eta = 0.1,
                            max_depth = 4,
                            subsample = 0.8,
                            colsample_bytree = 0.8,
                            min_child_weight = 1,
                            nrounds = 500,
                            early_stopping_rounds = 20,
                            training_size = 0.8,
                            eps = 1e-6,
                            seed = NULL,
                            return_fit = FALSE) {
  # check xgboost availability
  if (!requireNamespace("xgboost", quietly = TRUE)) stop("Install xgboost first.")
  # set the seed once so the point estimate and bootstrap are reproducible
  if (!is.null(seed)) set.seed(seed)
  # add equal reference weights if absent
  if (!(ref_weight_var %in% names(reference_survey))) reference_survey[[ref_weight_var]] <- 1
  # infer selection variables from shared columns
  if (is.null(selection_vars)) {
    selection_vars <- setdiff(intersect(names(data), names(reference_survey)),
                              c(outcome_var, ref_weight_var, "count"))
  }
  # stop if no selection variables are found
  if (length(selection_vars) == 0) stop("No shared selection variables found.")
  # check outcome variable
  if (!(outcome_var %in% names(data))) stop("Missing in data: ", outcome_var)
  # clean non-probability sample
  data0 <- data |>
    # keep complete rows
    filter(if_all(all_of(c(selection_vars, outcome_var)), ~ !is.na(.x))) |>
    # harmonise selection variables
    mutate(across(all_of(selection_vars), as.character))
  # clean reference sample
  ref0 <- reference_survey |>
    # keep complete rows
    filter(if_all(all_of(c(selection_vars, ref_weight_var)), ~ !is.na(.x))) |>
    # harmonise selection variables
    mutate(across(all_of(selection_vars), as.character))
  # drop unsupported sample levels
  for (v in selection_vars) {
    # keep only levels present in reference
    data0 <- data0[data0[[v]] %in% unique(ref0[[v]]), , drop = FALSE]
  }
  # stop if no data remains
  if (nrow(data0) == 0) stop("No data rows remain after matching reference support.")
  # reset row names before bootstrap indexing
  rownames(data0) <- NULL
  # one-hot encode predictors with aligned columns (reuse the NN helper)
  enc <- make_nn_dummy_matrices(data0, ref0, selection_vars)
  # non-prob design matrix
  x_np <- enc$x_data
  # reference design matrix
  x_ref <- enc$x_pop
  # reference weights aligned to x_ref rows
  ref_w <- enc$population[[ref_weight_var]]
  # make outcome categorical
  y <- factor(enc$data[[outcome_var]])
  # collect outcome levels
  outcome_levels <- levels(y)
  # helper: fit propensity model and return Hajek weighted shares
  est_once <- function(x_np, y, x_ref, ref_w) {
    # fit the gradient-boosted selection model on the pooled samples
    res <- fit_ipw_xgb_once(x_np, x_ref, ref_w,
                            eta = eta, max_depth = max_depth, subsample = subsample,
                            colsample_bytree = colsample_bytree, min_child_weight = min_child_weight,
                            nrounds = nrounds, early_stopping_rounds = early_stopping_rounds,
                            training_size = training_size)
    # clip propensities away from 0/1 for numerical stability
    p <- pmin(pmax(res$p, eps), 1 - eps)
    # inverse-odds pseudo-weights (population/non-prob density ratio)
    wt <- (1 - p) / p
    # self-normalised (Hajek) weighted shares per outcome category
    est <- vapply(outcome_levels, function(k) sum(wt * (y == k)) / sum(wt), numeric(1))
    # return a named estimate vector
    setNames(est, outcome_levels)
  }
  # compute the point estimate on the full samples
  estimate <- est_once(x_np, y, x_ref, ref_w)
  # bootstrap standard errors if requested
  se <- if (var_method == "bootstrap" && num_boot > 1) {
    # refit on resampled non-prob and reference rows
    boot <- replicate(num_boot, {
      # resample non-prob rows with replacement
      i <- sample(nrow(x_np), nrow(x_np), replace = TRUE)
      # resample reference rows with replacement
      j <- sample(nrow(x_ref), nrow(x_ref), replace = TRUE)
      # bootstrap estimate
      est_once(x_np[i, , drop = FALSE], y[i], x_ref[j, , drop = FALSE], ref_w[j])
    })
    # category-level standard errors across bootstrap replicates
    apply(boot, 1, sd, na.rm = TRUE)
  } else {
    # analytic variance is not available for the boosted model
    if (var_method != "bootstrap") message("Only var_method = 'bootstrap' is supported; returning NA standard errors.")
    # no standard errors
    rep(NA_real_, length(estimate))
  }
  # normal critical value
  z <- qnorm(1 - (1 - conf_level) / 2)
  # build clean estimate table (same shape as nonprob_ipw)
  out <- data.frame(
    # outcome category
    category = outcome_levels,
    # adjusted estimate
    estimate = as.numeric(estimate),
    # standard error
    se = as.numeric(se),
    # confidence interval lower bound
    ci_low = as.numeric(estimate - z * se),
    # confidence interval upper bound
    ci_high = as.numeric(estimate + z * se),
    # clean row names
    row.names = NULL
  )
  # return fit if requested
  if (return_fit) return(list(estimates = out, fit = est_once(x_np, y, x_ref, ref_w), data = data0, reference = ref0))
  # return estimates only
  out
}
