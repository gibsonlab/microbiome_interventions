#install.packages("devtools")
#install.packages("tidyverse")
#devtools::install_github("gathanei/xyz")
#devtools::install_github("krisrs1128/mbtransfer")

library(mbtransfer)
library(tidyverse)
library(glue)

set.seed(41)

rmse <- function(actual, predicted) {sqrt(mean((actual - predicted)^2, na.rm = TRUE))}

rmse_vec <- function(actual, predicted) {squared_errors <- (actual - predicted)^2
  mean_squared_errors <- rowMeans(squared_errors, na.rm = TRUE)
  sqrt(mean_squared_errors)}

subject <- read_csv("../_data/healthy/subject.csv")
interventions <- read_csv("../_data/healthy/interventions.csv") |> column_to_rownames("sample")
reads <- read_csv("../_data/healthy/reads.csv") |> column_to_rownames("sample")
samples <- read_csv("../_data/healthy/samples.csv")
qpcr <- read_tsv("../_data/healthy/qpcr.tsv") |> column_to_rownames("sample")


qpcr$median <- apply(qpcr, 1, median, na.rm=T)
qpcr_medians<- as.matrix(qpcr$median)

R <- as.matrix(reads)
R <- R+.1 # we choose .1 instead 1 so that minimum relative abundance of a zero read is always 1 order of magnitude lower than limit of detection in each sample
R <- R / rowSums(R) 
R <- t(t(R) * qpcr$median)
R <- log(R)

print(paste("max R=", max(R) ))
print(paste("min R=", min(R) ))

if (any(R<0)){warning("negative values being fed to model probably not a good idea with quadratic interactions")}

ts <- R |>
  ts_from_dfs(interventions, samples, subject) |>
  interpolate(method = "linear")

# create another time series that can easily identify 0 read samples so we make a raw reads time series for later use (not for training)
R_reads <- as.matrix(reads)
ts_reads <- R_reads |>
  ts_from_dfs(interventions, samples, subject) |>
  interpolate(method = "linear")

fits <- list()
ts_preds <- list()
subject_data(ts)

# mice are m2, m3, m4, m5 ... so  here we train m2,m3,m4 and then tests on m5, starting from only knowing the initial condition

ts_missing <- subset_values(ts, 1)

fits[["hold-out-m5"]] <- mbtransfer(ts[c(1:3)], P = 2, Q = 2)
ts_preds[["hold-out-m5"]] <- predict(fits[["hold-out-m5"]], ts_missing[c(4)])

# check forecasting error example 
y_predict<-ts_preds[["hold-out-m5"]]@series[["m5"]]@values
y_truth<-ts@series[["m5"]]@values
y_truth_reads<-ts_reads@series[["m5"]]@values

nr <- nrow(y_truth_reads)
y_truth_reads[1:nr, c(12:13, 15,17, 19:20, 26:27,34, 40:41, 48:49, 55:56)] <- 0

any(y_predict<0)
sum(y_predict<0)
if (any(y_predict<0)){warning("negative abundance values in forecast - but these will be exponentiated not so bad, just bad for model")}

# data is transformed so we want to exponentiate it and renormalize for a better comparison
y_predict_rel<- exp(y_predict) / colSums(exp(y_predict))
y_truth_rel<-exp(y_truth)/colSums(exp(y_truth))

any(y_predict_rel<0)
sum(y_predict_rel<0)
if (any(y_predict_rel<0)){warning("negative relative abundances - can not be evaluated directly")}

# compute error with log10 relative abundnces to compare all methods on equal footing (all trained on different scales)
rmse_out<-rmse(log10(y_predict_rel), log10(y_truth_rel))
print(rmse_out)

rmse_vec_out<-rmse_vec(log10(y_predict_rel), log10(y_truth_rel))
print(median(rmse_vec_out, na.rm = TRUE))

# only calculate error for the taxa-timepoints that have nonzero reads and remove interpolated timepoints like we do in the analysis in our paper
 
y_truth_rel_non_zero<-y_truth_rel
y_truth_rel_non_zero[y_truth_reads ==0] <- NA #mask values that originally had zero reads 

rmse_out<-rmse(log10(y_predict_rel), log10(y_truth_rel_non_zero))
print(rmse_out)

rmse_vec_out<-rmse_vec(log10(y_predict_rel), log10(y_truth_rel_non_zero))
print(median(rmse_vec_out, na.rm = TRUE))



################################################################################
#repeat the analysis but with slightly reduced scale but don't pass negative values !!!

R <- as.matrix(reads)
R <- R+.1 # we choose .1 instead 1 so that minimum relative abundance of a zero read is always 1 order of magnitude lower than limit of detection in each sample
R <- R / rowSums(R) 
R <- t(t(R) * qpcr$median)/1e3
R <- log(R)

print(paste("max R=", max(R) ))
print(paste("min R=", min(R) ))

if (any(R<0)){warning("negative values being fed to model probably not a good idea with quadratic interactions")}

ts <- R |>
  ts_from_dfs(interventions, samples, subject) |>
  interpolate(method = "linear")

fits <- list()
ts_preds <- list()
subject_data(ts)

# mice are m2, m3, m4, m5 ... so  here we train m2,m3,m4 and then tests on m5, starting from only knowing the initial condition

ts_missing <- subset_values(ts, 1)

fits[["hold-out-m5"]] <- mbtransfer(ts[c(1:3)], P = 2, Q = 2)
ts_preds[["hold-out-m5"]] <- predict(fits[["hold-out-m5"]], ts_missing[c(4)])

# check forecasting error example 
y_predict<-ts_preds[["hold-out-m5"]]@series[["m5"]]@values
y_truth<-ts@series[["m5"]]@values
y_truth_reads<-ts_reads@series[["m5"]]@values

nr <- nrow(y_truth_reads)
y_truth_reads[1:nr, c(12:13, 15,17, 19:20, 26:27,34, 40:41, 48:49, 55:56)] <- 0

any(y_predict<0)
sum(y_predict<0)
if (any(y_predict<0)){warning("negative abundance values in forecast - but these will be exponentiated not so bad, just bad for model")}

# data is transformed so we want to exponentiate it and renormalize for a better comparison
y_predict_rel<- exp(y_predict) / colSums(exp(y_predict))
y_truth_rel<-exp(y_truth)/colSums(exp(y_truth))

any(y_predict_rel<0)
sum(y_predict_rel<0)
if (any(y_predict_rel<0)){warning("negative relative abundances - can not be evaluated directly")}

# compute error with log10 relative abundnces to compare all methods on equal footing (all trained on different scales)
rmse_out<-rmse(log10(y_predict_rel), log10(y_truth_rel))
print(rmse_out)

rmse_vec_out<-rmse_vec(log10(y_predict_rel), log10(y_truth_rel))
print(median(rmse_vec_out, na.rm = TRUE))

# only calculate error for the taxa-timepoints that have nonzero reads and remove interpolated timepoints like we do in the analysis in our paper

y_truth_rel_non_zero<-y_truth_rel
y_truth_rel_non_zero[y_truth_reads ==0] <- NA #mask values that originally had zero reads 

rmse_out<-rmse(log10(y_predict_rel), log10(y_truth_rel_non_zero))
print(rmse_out)

rmse_vec_out<-rmse_vec(log10(y_predict_rel), log10(y_truth_rel_non_zero))
print(median(rmse_vec_out, na.rm = TRUE))

#modest improvemebnt but cant scale any more aggressivly otherwise we have negative values after log transform

##################
# what about training on raw CFU/g [crashes during inference - unless we reduce the stepsize?]

R <- as.matrix(reads)
R <- R+.1 # we choose .1 instead 1 so that minimum relative abundance of a zero read is always 1 order of magnitude lower than limit of detection in each sample
R <- R / rowSums(R) 
R <- t(t(R) * qpcr$median)
#R <- log(R) # data is in CFU/g no nonlinear transform

print(paste("max R=", max(R) ))
print(paste("min R=", min(R) ))

if (any(R<0)){warning("negative values being fed to model probably not a good idea with quadratic interactions")}


ts <- R |>
  ts_from_dfs(interventions, samples, subject) |>
  interpolate(method = "linear")

fits <- list()
ts_preds <- list()
subject_data(ts)

# mice are m2, m3, m4, m5 ... so  here we train m2,m3,m4 and then tests on m5, starting from only knowing the initial condition

ts_missing <- subset_values(ts, 1)

fits[["hold-out-m5"]] <- mbtransfer(ts[c(1:3)], P = 2, Q = 2) # model crashes during inference
ts_preds[["hold-out-m5"]] <- predict(fits[["hold-out-m5"]], ts_missing[c(4)])

# check forecasting error example 
y_predict<-ts_preds[["hold-out-m5"]]@series[["m5"]]@values
y_truth<-ts@series[["m5"]]@values
y_truth_reads<-ts_reads@series[["m5"]]@values

nr <- nrow(y_truth_reads)
y_truth_reads[1:nr, c(12:13, 15,17, 19:20, 26:27,34, 40:41, 48:49, 55:56)] <- 0

any(y_predict<0)
sum(y_predict<0)
if (any(y_predict<0)){warning("negative abundance values in forecast - these are in rescaled CFU/g so negative values are a problem")}

# set negative abundances to min(R)
y_predict[y_predict<0]<-min(R)/100
sum(y_predict<0)

# data is transformed so we want to exponentiate it and renormalize for a better comparison (log relasumtive abundance)
y_predict_rel<-(y_predict) #/ colSums((y_predict))
y_truth_rel<-(y_truth)#/colSums((y_truth))

rmse_out<-rmse(log10(y_predict_rel), log10(y_truth_rel))
print(rmse_out)

rmse_vec_out<-rmse_vec(log10(y_predict_rel), log10(y_truth_rel))
print(median(rmse_vec_out, na.rm = TRUE))

# only calculate error for the taxa-timepoints that have nonzero reads like we do in the analysis in our paper

y_truth_rel_non_zero<-y_truth_rel
y_truth_rel_non_zero[y_truth_reads == 0 ] <- NA #mask values that originally had zero reads 

rmse_out<-rmse(log10(y_predict_rel), log10(y_truth_rel_non_zero))
print(rmse_out)

rmse_vec_out<-rmse_vec(log10(y_predict_rel), log10(y_truth_rel_non_zero))
print(median(rmse_vec_out, na.rm = TRUE))

