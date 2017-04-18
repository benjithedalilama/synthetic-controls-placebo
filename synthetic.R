setwd('.')
source('synth_control_functions.R')

packages <- c('parallel', 'splines', 'zoo', 'lubridate', 'RcppRoll', 'CausalImpact', 'ggplot2', 'reshape')
packageHandler(packages, update_packages=TRUE, install_packages=TRUE)
sapply(packages, library, quietly = TRUE, character.only = TRUE)

#Detects number of available cores on computers. Used for parallel processing to speed up analysis.
n_cores <- detectCores()
set.seed(1)
prelog_data <- read.csv("Dataset S1 Brazil.csv")
group_name <- 'age_group'
groups <- as.character(unique(unlist(prelog_data[, group_name], use.names = FALSE)))
if (exists('exclude_group')) {groups <- groups[!(groups %in% exclude_group)]}

date_name <- 'date'
start_date <- as.Date("2003-01-01")
end_date <- as.Date("2013-12-01")
prelog_data[, date_name] <- formatDate(prelog_data[, date_name])
outcome_name <- 'J12_18'
denom_name <- 'A10_B99_nopneumo'
prelog_data <- setNames(lapply(groups, FUN = splitGroup, ungrouped_data = prelog_data, group_name = group_name, date_name = date_name, start_date = start_date, end_date = end_date, no_filter = c(group_name, date_name, outcome_name, denom_name)), groups)
#if (exists('exclude_group')) {prelog_data <- prelog_data[!(names(prelog_data) %in% exclude_group)]}

#Log-transform all variables, adding 0.5 to counts of 0.
ds <- setNames(lapply(prelog_data, FUN = logTransform, no_log = c(group_name, date_name)), groups)
time_points <- unique(ds[[1]][, date_name])

ds <- lapply(ds, function(ds) {
	if (!(denom_name %in% colnames(ds))) {
		ds[denom_name] <- 0
	}
	return(ds)
})

exclude_covar <- c()
sparse_groups <- sapply(ds, function(ds) {
	return(ncol(ds[!(colnames(ds) %in% c(date_name, group_name, denom_name, outcome_name, exclude_covar))]) == 0)
})
ds <- ds[!sparse_groups]
groups <- groups[!sparse_groups]

code_change <- TRUE

intervention_date <- as.Date('2008-01-01')
#Process and standardize the covariates. For the Brazil data, adjust for 2008 coding change.
covars_full <- setNames(lapply(ds, makeCovars, code_change = code_change, intervention_date = intervention_date, time_points = time_points), groups)
covars_full <- sapply(covars_full, FUN = function(covars) {covars[, !(colnames(covars) %in% exclude_covar), drop = FALSE]})
covars_time <- setNames(lapply(covars_full, FUN = function(covars) {as.data.frame(list(time_index = 1:nrow(covars)))}), groups)

#Standardize the outcome variable and save the original mean and SD for later analysis.
outcome      <- sapply(ds, FUN = function(data) {scale(data[, outcome_name])})
outcome_mean <- sapply(ds, FUN = function(data) {mean(data[, outcome_name])})
outcome_sd   <- sapply(ds, FUN = function(data) {sd(data[, outcome_name])})
outcome_plot <- exp(t(t(outcome) * outcome_sd + outcome_mean))
outcome_offset      <- sapply(ds, FUN = function(data) {data[, outcome_name] - data[, denom_name]})
outcome_offset_mean <- colMeans(outcome_offset)
outcome_offset_sd   <- sapply(ds, FUN = function(data) {sd(data[, outcome_name] - data[, denom_name])})
outcome_offset      <- scale(outcome_offset)

#Combine the outcome, covariates, and time point information.
data_full <- setNames(lapply(groups, makeTimeSeries, outcome = outcome,        covars = covars_full, time_points = time_points), groups)
data_time <- setNames(lapply(groups, makeTimeSeries, outcome = outcome_offset, covars = covars_time, time_points = time_points), groups)

n_cores<-4
n_seasons <- 12
cl <- makeCluster(n_cores)
clusterEvalQ(cl, {library(CausalImpact, quietly = TRUE); library(lubridate, quietly = TRUE)})
clusterExport(cl, c('doCausalImpact', 'impactExtract', 'intervention_date', 'time_points', 'n_seasons'), environment())
impact_full <- setNames(parLapply(cl, data_full, doCausalImpact, intervention_date = intervention_date, time_points = time_points, n_seasons = n_seasons), groups)
#impact_time <- setNames(parLapply(cl, data_time, doCausalImpact, intervention_date = intervention_date, time_points = time_points, n_seasons = n_seasons, trend = TRUE), groups)

stopCluster(cl)

plot.new()
plot(impact_full$`0`) + ggtitle(paste(intervention_date, 'Brazil age group 0'))
plot(impact_full$`1`) + ggtitle(paste(intervention_date, 'Brazil age group 1'))
plot(impact_full$`2`) + ggtitle(paste(intervention_date, 'Brazil age group 2'))
plot(impact_full$`3`) + ggtitle(paste(intervention_date, 'Brazil age group 3'))
plot(impact_full$`4`) + ggtitle(paste(intervention_date, 'Brazil age group 4'))
plot(impact_full$`5`) + ggtitle(paste(intervention_date, 'Brazil age group 5'))
plot(impact_full$`6`) + ggtitle(paste(intervention_date, 'Brazil age group 6'))
plot(impact_full$`7`) + ggtitle(paste(intervention_date, 'Brazil age group 7'))
plot(impact_full$`8`) + ggtitle(paste(intervention_date, 'Brazil age group 8'))
plot(impact_full$`9`) + ggtitle(paste(intervention_date, 'Brazil age group 9'))
paste(impact_full$`0`$summary$RelEffect,'[',impact_full$`0`$summary$RelEffect.lower,',',impact_full$`0`$summary$RelEffect.upper,']',impact_full$`0`$summary$p)
paste(impact_full$`1`$summary$RelEffect,'[',impact_full$`1`$summary$RelEffect.lower,',',impact_full$`1`$summary$RelEffect.upper,']',impact_full$`1`$summary$p)
paste(impact_full$`2`$summary$RelEffect,'[',impact_full$`2`$summary$RelEffect.lower,',',impact_full$`2`$summary$RelEffect.upper,']',impact_full$`2`$summary$p)
paste(impact_full$`3`$summary$RelEffect,'[',impact_full$`3`$summary$RelEffect.lower,',',impact_full$`3`$summary$RelEffect.upper,']',impact_full$`3`$summary$p)
paste(impact_full$`4`$summary$RelEffect,'[',impact_full$`4`$summary$RelEffect.lower,',',impact_full$`4`$summary$RelEffect.upper,']',impact_full$`4`$summary$p)
paste(impact_full$`5`$summary$RelEffect,'[',impact_full$`5`$summary$RelEffect.lower,',',impact_full$`5`$summary$RelEffect.upper,']',impact_full$`5`$summary$p)
paste(impact_full$`6`$summary$RelEffect,'[',impact_full$`6`$summary$RelEffect.lower,',',impact_full$`6`$summary$RelEffect.upper,']',impact_full$`6`$summary$p)
paste(impact_full$`7`$summary$RelEffect,'[',impact_full$`7`$summary$RelEffect.lower,',',impact_full$`7`$summary$RelEffect.upper,']',impact_full$`7`$summary$p)
paste(impact_full$`8`$summary$RelEffect,'[',impact_full$`8`$summary$RelEffect.lower,',',impact_full$`8`$summary$RelEffect.upper,']',impact_full$`8`$summary$p)
paste(impact_full$`9`$summary$RelEffect,'[',impact_full$`9`$summary$RelEffect.lower,',',impact_full$`9`$summary$RelEffect.upper,']',impact_full$`9`$summary$p)

intervention_date <- as.Date('2009-01-01')
#Process and standardize the covariates. For the Brazil data, adjust for 2008 coding change.
covars_full <- setNames(lapply(ds, makeCovars, code_change = code_change, intervention_date = intervention_date, time_points = time_points), groups)
covars_full <- sapply(covars_full, FUN = function(covars) {covars[, !(colnames(covars) %in% exclude_covar), drop = FALSE]})
covars_time <- setNames(lapply(covars_full, FUN = function(covars) {as.data.frame(list(time_index = 1:nrow(covars)))}), groups)

#Standardize the outcome variable and save the original mean and SD for later analysis.
outcome      <- sapply(ds, FUN = function(data) {scale(data[, outcome_name])})
outcome_mean <- sapply(ds, FUN = function(data) {mean(data[, outcome_name])})
outcome_sd   <- sapply(ds, FUN = function(data) {sd(data[, outcome_name])})
outcome_plot <- exp(t(t(outcome) * outcome_sd + outcome_mean))
outcome_offset      <- sapply(ds, FUN = function(data) {data[, outcome_name] - data[, denom_name]})
outcome_offset_mean <- colMeans(outcome_offset)
outcome_offset_sd   <- sapply(ds, FUN = function(data) {sd(data[, outcome_name] - data[, denom_name])})
outcome_offset      <- scale(outcome_offset)

#Combine the outcome, covariates, and time point information.
data_full <- setNames(lapply(groups, makeTimeSeries, outcome = outcome,        covars = covars_full, time_points = time_points), groups)
data_time <- setNames(lapply(groups, makeTimeSeries, outcome = outcome_offset, covars = covars_time, time_points = time_points), groups)

n_cores<-4
n_seasons <- 12
cl <- makeCluster(n_cores)
clusterEvalQ(cl, {library(CausalImpact, quietly = TRUE); library(lubridate, quietly = TRUE)})
clusterExport(cl, c('doCausalImpact', 'impactExtract', 'intervention_date', 'time_points', 'n_seasons'), environment())
impact_full <- setNames(parLapply(cl, data_full, doCausalImpact, intervention_date = intervention_date, time_points = time_points, n_seasons = n_seasons), groups)
#impact_time <- setNames(parLapply(cl, data_time, doCausalImpact, intervention_date = intervention_date, time_points = time_points, n_seasons = n_seasons, trend = TRUE), groups)

stopCluster(cl)

plot.new()
plot(impact_full$`0`) + ggtitle(paste(intervention_date, 'Brazil age group 0'))
plot(impact_full$`1`) + ggtitle(paste(intervention_date, 'Brazil age group 1'))
plot(impact_full$`2`) + ggtitle(paste(intervention_date, 'Brazil age group 2'))
plot(impact_full$`3`) + ggtitle(paste(intervention_date, 'Brazil age group 3'))
plot(impact_full$`4`) + ggtitle(paste(intervention_date, 'Brazil age group 4'))
plot(impact_full$`5`) + ggtitle(paste(intervention_date, 'Brazil age group 5'))
plot(impact_full$`6`) + ggtitle(paste(intervention_date, 'Brazil age group 6'))
plot(impact_full$`7`) + ggtitle(paste(intervention_date, 'Brazil age group 7'))
plot(impact_full$`8`) + ggtitle(paste(intervention_date, 'Brazil age group 8'))
plot(impact_full$`9`) + ggtitle(paste(intervention_date, 'Brazil age group 9'))
paste(impact_full$`0`$summary$RelEffect,'[',impact_full$`0`$summary$RelEffect.lower,',',impact_full$`0`$summary$RelEffect.upper,']',impact_full$`0`$summary$p)
paste(impact_full$`1`$summary$RelEffect,'[',impact_full$`1`$summary$RelEffect.lower,',',impact_full$`1`$summary$RelEffect.upper,']',impact_full$`1`$summary$p)
paste(impact_full$`2`$summary$RelEffect,'[',impact_full$`2`$summary$RelEffect.lower,',',impact_full$`2`$summary$RelEffect.upper,']',impact_full$`2`$summary$p)
paste(impact_full$`3`$summary$RelEffect,'[',impact_full$`3`$summary$RelEffect.lower,',',impact_full$`3`$summary$RelEffect.upper,']',impact_full$`3`$summary$p)
paste(impact_full$`4`$summary$RelEffect,'[',impact_full$`4`$summary$RelEffect.lower,',',impact_full$`4`$summary$RelEffect.upper,']',impact_full$`4`$summary$p)
paste(impact_full$`5`$summary$RelEffect,'[',impact_full$`5`$summary$RelEffect.lower,',',impact_full$`5`$summary$RelEffect.upper,']',impact_full$`5`$summary$p)
paste(impact_full$`6`$summary$RelEffect,'[',impact_full$`6`$summary$RelEffect.lower,',',impact_full$`6`$summary$RelEffect.upper,']',impact_full$`6`$summary$p)
paste(impact_full$`7`$summary$RelEffect,'[',impact_full$`7`$summary$RelEffect.lower,',',impact_full$`7`$summary$RelEffect.upper,']',impact_full$`7`$summary$p)
paste(impact_full$`8`$summary$RelEffect,'[',impact_full$`8`$summary$RelEffect.lower,',',impact_full$`8`$summary$RelEffect.upper,']',impact_full$`8`$summary$p)
paste(impact_full$`9`$summary$RelEffect,'[',impact_full$`9`$summary$RelEffect.lower,',',impact_full$`9`$summary$RelEffect.upper,']',impact_full$`9`$summary$p)

intervention_date <- as.Date('2010-01-01')
#Process and standardize the covariates. For the Brazil data, adjust for 2008 coding change.
covars_full <- setNames(lapply(ds, makeCovars, code_change = code_change, intervention_date = intervention_date, time_points = time_points), groups)
covars_full <- sapply(covars_full, FUN = function(covars) {covars[, !(colnames(covars) %in% exclude_covar), drop = FALSE]})
covars_time <- setNames(lapply(covars_full, FUN = function(covars) {as.data.frame(list(time_index = 1:nrow(covars)))}), groups)

#Standardize the outcome variable and save the original mean and SD for later analysis.
outcome      <- sapply(ds, FUN = function(data) {scale(data[, outcome_name])})
outcome_mean <- sapply(ds, FUN = function(data) {mean(data[, outcome_name])})
outcome_sd   <- sapply(ds, FUN = function(data) {sd(data[, outcome_name])})
outcome_plot <- exp(t(t(outcome) * outcome_sd + outcome_mean))
outcome_offset      <- sapply(ds, FUN = function(data) {data[, outcome_name] - data[, denom_name]})
outcome_offset_mean <- colMeans(outcome_offset)
outcome_offset_sd   <- sapply(ds, FUN = function(data) {sd(data[, outcome_name] - data[, denom_name])})
outcome_offset      <- scale(outcome_offset)

#Combine the outcome, covariates, and time point information.
data_full <- setNames(lapply(groups, makeTimeSeries, outcome = outcome,        covars = covars_full, time_points = time_points), groups)
data_time <- setNames(lapply(groups, makeTimeSeries, outcome = outcome_offset, covars = covars_time, time_points = time_points), groups)

n_cores<-4
n_seasons <- 12
cl <- makeCluster(n_cores)
clusterEvalQ(cl, {library(CausalImpact, quietly = TRUE); library(lubridate, quietly = TRUE)})
clusterExport(cl, c('doCausalImpact', 'impactExtract', 'intervention_date', 'time_points', 'n_seasons'), environment())
impact_full <- setNames(parLapply(cl, data_full, doCausalImpact, intervention_date = intervention_date, time_points = time_points, n_seasons = n_seasons), groups)
#impact_time <- setNames(parLapply(cl, data_time, doCausalImpact, intervention_date = intervention_date, time_points = time_points, n_seasons = n_seasons, trend = TRUE), groups)

stopCluster(cl)

plot.new()
plot(impact_full$`0`) + ggtitle(paste(intervention_date, 'Brazil age group 0'))
plot(impact_full$`1`) + ggtitle(paste(intervention_date, 'Brazil age group 1'))
plot(impact_full$`2`) + ggtitle(paste(intervention_date, 'Brazil age group 2'))
plot(impact_full$`3`) + ggtitle(paste(intervention_date, 'Brazil age group 3'))
plot(impact_full$`4`) + ggtitle(paste(intervention_date, 'Brazil age group 4'))
plot(impact_full$`5`) + ggtitle(paste(intervention_date, 'Brazil age group 5'))
plot(impact_full$`6`) + ggtitle(paste(intervention_date, 'Brazil age group 6'))
plot(impact_full$`7`) + ggtitle(paste(intervention_date, 'Brazil age group 7'))
plot(impact_full$`8`) + ggtitle(paste(intervention_date, 'Brazil age group 8'))
plot(impact_full$`9`) + ggtitle(paste(intervention_date, 'Brazil age group 9'))
paste(impact_full$`0`$summary$RelEffect,'[',impact_full$`0`$summary$RelEffect.lower,',',impact_full$`0`$summary$RelEffect.upper,']',impact_full$`0`$summary$p)
paste(impact_full$`1`$summary$RelEffect,'[',impact_full$`1`$summary$RelEffect.lower,',',impact_full$`1`$summary$RelEffect.upper,']',impact_full$`1`$summary$p)
paste(impact_full$`2`$summary$RelEffect,'[',impact_full$`2`$summary$RelEffect.lower,',',impact_full$`2`$summary$RelEffect.upper,']',impact_full$`2`$summary$p)
paste(impact_full$`3`$summary$RelEffect,'[',impact_full$`3`$summary$RelEffect.lower,',',impact_full$`3`$summary$RelEffect.upper,']',impact_full$`3`$summary$p)
paste(impact_full$`4`$summary$RelEffect,'[',impact_full$`4`$summary$RelEffect.lower,',',impact_full$`4`$summary$RelEffect.upper,']',impact_full$`4`$summary$p)
paste(impact_full$`5`$summary$RelEffect,'[',impact_full$`5`$summary$RelEffect.lower,',',impact_full$`5`$summary$RelEffect.upper,']',impact_full$`5`$summary$p)
paste(impact_full$`6`$summary$RelEffect,'[',impact_full$`6`$summary$RelEffect.lower,',',impact_full$`6`$summary$RelEffect.upper,']',impact_full$`6`$summary$p)
paste(impact_full$`7`$summary$RelEffect,'[',impact_full$`7`$summary$RelEffect.lower,',',impact_full$`7`$summary$RelEffect.upper,']',impact_full$`7`$summary$p)
paste(impact_full$`8`$summary$RelEffect,'[',impact_full$`8`$summary$RelEffect.lower,',',impact_full$`8`$summary$RelEffect.upper,']',impact_full$`8`$summary$p)
paste(impact_full$`9`$summary$RelEffect,'[',impact_full$`9`$summary$RelEffect.lower,',',impact_full$`9`$summary$RelEffect.upper,']',impact_full$`9`$summary$p)