#'Predict a biomarker based on the other using specified model
#'
#' @param dat The data frame after transformation
#' @param studies A column for all study IDs
#' @param mod_obj A model object
#' @param N Number of simulations
#' @param seed Seed value for reproducibility; use default if not provided
#' @return A vector of predicted probabilities
#' @examples
#' predicted_a1c <- cw_predict(data, data$study_id, mod_selected, N = 2000)
cw_predict <- function(dat, studies, mod_obj, N = 2000, seed = 2022) {
    # Fixed effect
    if (N > 5000) stop('Cannot run more than 5000 iterations')
    fixef <- mod_obj$fixefs[1:N,]
    boot <- fixef %*% t(as.matrix(dat))
    unique_studies <- unique(studies)

    set.seed(seed)

    # Study RE - individuals from the same study get the same REs
    for (i in 1:length(unique_studies)) {
        ranef_i <- rnorm(n = N, mean = 0, sd = mod_obj$ranef_sd)
        boot[,studies == unique_studies[i]] <- boot[,studies == unique_studies[i]] + ranef_i
    }

    boot <- t(arm::invlogit(boot))

    return(boot)
}

#'Predict missing HbA1c based on FPG using an appropriate model
#'
#' @param d The data frame after transformation
#' @param all_model_coefs An list of all models
#' @param N_sim Number of simulations
#' @param seed Seed value for reproducibility; use default if not provided
#' @param return_draws logical value: if draw level prediction is returned as a matrix
#' @return if return_draws is TRUE, a list including a vector of predicted probabilities
#'         and a matrix of predictions at draw level; if FALSE, only a vector
#' @examples
#' predictions <- transform(data, models, N_sim = 2000)
#' @export
predict_missing_hba1c <- function(d, all_model_coefs, N_sim = 2000, seed = 2022, return_draws = FALSE, verbose = TRUE) {

    # create indicators for cross-walk
    to_predict <- d$medication == 0
    to_predict_a1c <- to_predict & is.na(d$a1c.var) & !is.na(d$fpg.var)

    has_reg <- d$Superregion != 'Oceania' & !is.na(d$Superregion)
    has_bmi <- !is.na(d$bmi.var) & is.finite(d$bmi.var)
    has_dev_fpg <- !is.na(d$handheld_fpg_cat)

    # Model groups based on availability of covariates
    group_fpgtoa1c <- factor(paste(as.integer(has_reg), as.integer(has_bmi), as.integer(has_dev_fpg)), levels = c('0 0 0', '0 0 1', '0 1 0', '0 1 1', '1 0 0', '1 0 1', '1 1 0', '1 1 1'))
    # Assign group names that match the model object
    # levels should be : 0 0 0   0 0 1   0 1 0   0 1 1   1 0 0   1 0 1   1 1 0   1 1 1
    model_groups <- c('noreg_none','noreg_nobmi','noreg_nodev','noreg_full','none','nobmi','nodev','full')
    levels(group_fpgtoa1c) <- model_groups

    # One-hot encoding for categorical variables
    encoding <- dataPreparation::build_encoding(d, cols = c("Superregion", "handheld_fpg_cat", "sex"), verbose = FALSE)
    dtmp <- dataPreparation::one_hot_encoder(d, encoding, verbose = FALSE)

    # Matching data columns with model coefficients
    dtmp <- dtmp %>%
        dplyr::select(study_id, sex.male, age.var, fpg.var, a1c.var, paste0('Superregion.', gsub('-|,','.',srs)), bmi.var, handheld.fpg.cat.Handheld) %>%
        rename(handheld_fpg_catHandheld = handheld.fpg.cat.Handheld,
               sexmale = sex.male)
    names(dtmp)[grep('Superregion', names(dtmp))] <- gsub('\\.', ',', gsub('Sub\\.', 'Sub-', gsub('High\\.','High-', gsub('Superregion\\.', 'Superregion', names(dtmp)[grep('Superregion',names(dtmp))]))))
    dtmp <- cbind(1, dtmp)
    names(dtmp)[1] <- '(Intercept)'

    # Adding interaction terms
    # fpg
    fpg_terms <- dtmp$fpg.var * dtmp %>% dplyr::select(starts_with('Superregion'))
    names(fpg_terms) <- paste0('fpg.var:', names(fpg_terms))
    d_fpg <- cbind(dtmp %>% dplyr::select(-contains('a1c')), fpg_terms)

    mat_predicted_a1c <- matrix(NA, nrow(data), N_sim)
    predicted_a1c <- rep(NA, nrow(d))

    # Predicting by model group
    for (grp in model_groups) {
        list1 <- which(to_predict_a1c & group_fpgtoa1c == grp)
        if (verbose) print(paste(grp, length(list1)))

        d_grp <- d_fpg[list1,]
        dx <- d[list1, ]

        if (length(list1)>0) {
            mod <- all_model_coefs$fpgtoa1c_mods[[paste0('mod_',grp)]]
            tmp <- d_grp %>% dplyr::select(all_of(colnames(mod$fixefs)))
            res <- cw_predict(tmp, d_grp$study_id, mod, N_sim)

            predicted_a1c[list1] <- matrixStats::rowMedians(res)
            mat_predicted_a1c[list1,] <- res
        }

    }

    if (return_draws) {
        return(list(
            predicted_a1c = predicted_a1c,
            mat_predicted_a1c = mat_predicted_a1c
        ))
    } else {
        return(predicted_a1c)
    }

}


#'Predict missing FPG based on HbA1c using an appropriate model
#'
#' @param d The data frame after transformation
#' @param all_model_coefs An list of all models
#' @param N_sim Number of simulations
#' @param seed Seed value for reproducibility; use default if not provided
#' @param return_draws logical value: if draw level prediction is returned as a matrix
#' @return if return_draws is TRUE, a list including a vector of predicted probabilities
#'         and a matrix of predictions at draw level; if FALSE, only a vector
#' @examples
#' predictions <- transform(data, models, N_sim = 2000)
#' @export
predict_missing_fpg <- function(d, all_model_coefs, N_sim = 2000, seed = 2022, return_draws = FALSE, verbose = TRUE) {

    # create indicators for cross-walk
    to_predict <- d$medication == 0
    to_predict_fpg <- to_predict & is.na(d$fpg.var) & !is.na(d$a1c.var)

    has_reg <- d$Superregion != 'Oceania' & !is.na(d$Superregion)
    has_bmi <- !is.na(d$bmi.var) & is.finite(d$bmi.var)
    has_dev_a1c <- !is.na(d$handheld_a1c_cat)

    # Model groups based on availability of covariates
    group_a1ctofpg <- factor(paste(as.integer(has_reg), as.integer(has_bmi), as.integer(has_dev_a1c)), levels = c('0 0 0', '0 0 1', '0 1 0', '0 1 1', '1 0 0', '1 0 1', '1 1 0', '1 1 1'))
    # Assign group names that match the model object
    # levels should be : 0 0 0   0 0 1   0 1 0   0 1 1   1 0 0   1 0 1   1 1 0   1 1 1
    model_groups <- c('noreg_none','noreg_nobmi','noreg_nodev','noreg_full','none','nobmi','nodev','full')
    levels(group_a1ctofpg) <- model_groups

    # One-hot encoding for categorical variables
    encoding <- dataPreparation::build_encoding(d, cols = c("Superregion", "handheld_a1c_cat", "sex"), verbose = FALSE)
    dtmp <- dataPreparation::one_hot_encoder(d, encoding, verbose = FALSE)

    # Matching data columns with model coefficients
    dtmp <- dtmp %>%
        dplyr::select(study_id, sex.male, age.var, fpg.var, a1c.var, paste0('Superregion.', gsub('-|,','.',srs)), bmi.var, handheld.a1c.cat.Handheld) %>%
        rename(handheld_a1c_catHandheld = handheld.a1c.cat.Handheld,
               sexmale = sex.male)
    names(dtmp)[grep('Superregion', names(dtmp))] <- gsub('\\.', ',', gsub('Sub\\.', 'Sub-', gsub('High\\.','High-', gsub('Superregion\\.', 'Superregion', names(dtmp)[grep('Superregion',names(dtmp))]))))
    dtmp <- cbind(1, dtmp)
    names(dtmp)[1] <- '(Intercept)'

    # Adding interaction terms
    # a1c
    a1c_terms <- dtmp$a1c.var * dtmp %>% dplyr::select(starts_with('Superregion'))
    names(a1c_terms) <- paste0('a1c.var:', names(a1c_terms))
    d_a1c <- cbind(dtmp %>% dplyr::select(-contains('fpg')), a1c_terms)

    mat_predicted_fpg <- matrix(NA, nrow(data), N_sim)
    predicted_fpg <- rep(NA, nrow(d))

    # Predicting by model group
    for (grp in model_groups) {
        list2 <- which(to_predict_fpg & group_a1ctofpg == grp)
        if (verbose) print(paste(grp, length(list2)))

        d_grp <- d_a1c[list2,]
        dx <- d[list2, ]

        if (length(list2)>0) {
            mod <- all_model_coefs$a1ctofpg_mods[[paste0('mod_',grp)]]
            tmp <- d_grp %>% dplyr::select(all_of(colnames(mod$fixefs)))
            res <- cw_predict(tmp, d_grp$study_id, mod, N_sim)

            predicted_fpg[list2] <- matrixStats::rowMedians(res)
            mat_predicted_fpg[list2,] <- res
        }

    }

    if (return_draws) {
        return(list(
            predicted_fpg = predicted_fpg,
            mat_predicted_fpg = mat_predicted_fpg
        ))
    } else  {
        return(predicted_fpg)
    }

}
