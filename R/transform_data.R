#'Transform data consistent with how models were fitted
#'
#' @param d The data frame with all necessary columns
#' @return A data frame with columns of standard names
#' @examples
#' data_transformed <- transform(data)
#' @export
transform_data <- function(d) {
    # check if columns all exist
    var_list <- c('iso', 'age', 'sex', 'study_id', 'medication')
    if (any(!var_list %in% names(d))) stop(paste(paste(var_list[!var_list%in%names(d)], collapse = ' '), 'not in the data frame'))
    pred_list <- c('fpg_clean', 'hba1c_clean')
    if (all(!pred_list %in% names(d))) stop(paste('None of', paste(pred_list, collapse = ' '), 'is in the data frame'))

    # Assign Superregion values if not already available
    if ((!'Superregion' %in% names(d)) | any(is.na(d$Superregion))) {
        # countrylist <- read.csv('country-list-2023-new.csv') # supplied within data.Rda
        if (!'Superregion' %in% names(d)) d$Superregion <- NA
        d$Superregion[is.na(d$Superregion)] <- countrylist$Superregion[match(d$iso[is.na(d$Superregion)], countrylist$iso)]
    }
    if (any(is.na(d$Superregion))) {
        stop(paste('ISO not recognised: correct ISO or provide Superregion directly\n', paste(unique(d$iso[is.na(d$Superregion)]), collapse = ', ')))
    }
    # historical renaming
    d$Superregion[which(d$Superregion == 'East and southeast Asia and the Pacific')] <- 'Southeast Asia, East Asia and the Pacific'

    # if Superregion was supplied but not in the standard names
    if (any(!d$Superregion %in% srs)) stop(paste('Super-regions:', paste(unique(d$Superregion[!d$Superregion %in% srs]), collapse = ', '), 'not in the standard super-region list'))
    d$Superregion <- factor(d$Superregion, levels = srs)

    # normalise variables using the same mean and SD as used in model development
    d$age.var <- (d$age - 50) / 15
    d$sex <- factor(d$sex, levels = c('female', 'male'))
    if ('hba1c_clean' %in% names(d)) d$a1c.var <- (d$hba1c_clean - 5.5) / 0.7 else d$a1c.var <- NA
    if ('fpg_clean' %in% names(d)) d$fpg.var <- (d$fpg_clean - 5.5) / 1 else d$fpg.var <- NA
    if ('bmi_clean' %in% names(d)) d$bmi.var <- (d$bmi_clean - 26.5) / 5 else d$bmi.var <- NA

    # recode device info dummy
    if (!'handheld_a1c' %in% names(d)) d$handheld_a1c <- NA
    if (!'handheld_fpg' %in% names(d)) d$handheld_fpg <- NA
    d$handheld_a1c_cat <- ifelse(is.na(d$handheld_a1c), NA, ifelse(d$handheld_a1c == 1, 'Handheld', 'Lab'))
    d$handheld_a1c_cat <- factor(d$handheld_a1c_cat, levels = c('Lab', 'Handheld'))
    d$handheld_fpg_cat <- ifelse(is.na(d$handheld_fpg), NA, ifelse(d$handheld_fpg == 1, 'Handheld', 'Lab'))
    d$handheld_fpg_cat <- factor(d$handheld_fpg_cat, levels = c('Lab', 'Handheld'))

    return(d[,c('sex','age.var','fpg.var','a1c.var','Superregion','handheld_fpg_cat','handheld_a1c_cat','study_id','bmi.var','medication')])
}
