#install.packages("devtools")
#install.packages("tidyverse")
#devtools::install_github("gathanei/xyz")
#devtools::install_github("krisrs1128/mbtransfer")


library(tidyverse)
library(glue)
library(mbtransfer)

eps<-1e-8

rmse <- function(actual, predicted) {sqrt(mean((actual - predicted)^2, na.rm = TRUE))}

rmse_vec <- function(actual, predicted) {
  squared_errors <- (actual - predicted)^2
  mean_squared_errors <- rowMeans(squared_errors, na.rm = TRUE)
  sqrt(mean_squared_errors)}

subject <- read_csv("../_data/healthy/subject.csv")
interventions <- read_csv("../_data/healthy/interventions.csv") |> column_to_rownames("sample")
reads <- read_csv("../_data/healthy/reads.csv") |> column_to_rownames("sample")
samples <- read_csv("../_data/healthy/samples.csv")

R <- as.matrix(reads)
R <- R+.1 # we choose .1 instead 1 so that minimum relative abundance of a zero read is always 1 order of magnitude lower than limit of detection in each sample
R <- R / rowSums(R) #relative abundance
R <- (R+min(R)/10) /min(R) # like assuming constant read depth for all samples with an epsilon value added back in 
# so we don't have 0 after log transform, this is equivalent to R/min(R)+.1 so its like normalized read depth
min(R) # this should 1.1, so that after log we still have nonzero abundance.

R <- log(R) # log transform there should be no negative values
min(R)
max(R)

if (any(R<0)){warning('negative values being fed to model probably not a good idea')}

ts <- R |>
  ts_from_dfs(interventions, samples, subject) |>
  interpolate(method = "linear")

# create another time series that we can easily identify 0 read samples
R_reads<- as.matrix(reads)
ts_reads <- R_reads |>
  ts_from_dfs(interventions, samples, subject) |>
  interpolate(method = "linear")

fits <- list()
ts_preds <- list()
subject_data(ts)

# mice are m2, m3, m4, m5 ... so  here we m2,m3,m4 and then tests on m5, starting from only knowing the initial condition

ts_missing <- subset_values(ts, 1)

fits[["hold-out-m5"]] <- mbtransfer(ts[c(1:3)], P = 2, Q = 2)
ts_preds[["hold-out-m5"]] <- predict(fits[["hold-out-m5"]], ts_missing[c(4)]) ## models trains but forecast crashes

y_predict<-ts_preds[["hold-out-m5"]]@series[["m5"]]@values
y_truth<-ts@series[["m5"]]@values
y_truth_reads<-ts_reads@series[["m5"]]@values

nr <- nrow(y_truth_reads)
y_truth_reads[1:nr, c(12:13, 15, 17, 19:20, 26:27, 34, 40:41, 48:49, 55:56)] <- 0 # interpolated columns ro be removed from analysis here - in the handoff youn deals with this directly

y_predict_rel <- exp(y_predict) / colSums(exp(y_predict))
y_truth_rel<-exp(y_truth)/colSums(exp(y_truth))

any(y_predict_rel<0)

out<-rmse_vec(log10(y_predict_rel+eps), log10(y_truth_rel+eps))
print(median(out, na.rm = TRUE))

out<-rmse(log10(y_predict_rel+eps), log10(y_truth_rel+eps))
print(out)

# only calculate error for the taxa-timepoints that have nonzero reads like we do in the analysis in our paper and remove interpolateed

y_truth_rel[y_truth_reads==0] <-NA

out<-rmse_vec(log10(y_predict_rel+eps), log10(y_truth_rel+eps))
print(median(out, na.rm = TRUE))

out<-rmse(log10(y_predict_rel+eps), log10(y_truth_rel+eps))
print(out)

# repeat analysis for saving to csv and handing off to Youn

## 5

fits[["hold-out-m5"]] <- mbtransfer(ts[c(1:3)], P = 2, Q = 2)
ts_preds[["hold-out-m5"]] <- predict(fits[["hold-out-m5"]], ts_missing[c(4)]) ## models trains but forecast crashes

y_predict<-ts_preds[["hold-out-m5"]]@series[["m5"]]@values
y_truth<-ts@series[["m5"]]@values
y_truth_reads<-ts_reads@series[["m5"]]@values

y_predict_rel <- exp(y_predict) / colSums(exp(y_predict))
y_truth_rel<-exp(y_truth)/colSums(exp(y_truth))
y_truth_rel[y_truth_reads==0] <- 0

y_predict_rel_5<-y_predict_rel
y_truth_rel_5<-y_truth_rel

## 4

fits[["hold-out-m4"]] <- mbtransfer(ts[c(1:2, 4)], P = 2, Q = 2)
ts_preds[["hold-out-m4"]] <- predict(fits[["hold-out-m4"]], ts_missing[c(3)])

y_predict<-ts_preds[["hold-out-m4"]]@series[["m4"]]@values
y_truth<-ts@series[["m4"]]@values
y_truth_reads<-ts_reads@series[["m4"]]@values

y_predict_rel <- exp(y_predict) / colSums(exp(y_predict))
y_truth_rel<-exp(y_truth)/colSums(exp(y_truth))
y_truth_rel[y_truth_reads==0] <- 0

y_predict_rel_4<-y_predict_rel
y_truth_rel_4<-y_truth_rel

## 3

fits[["hold-out-m3"]] <- mbtransfer(ts[c(1, 3:4)], P = 2, Q = 2)
ts_preds[["hold-out-m3"]] <- predict(fits[["hold-out-m3"]], ts_missing[c(2)])

y_predict<-ts_preds[["hold-out-m3"]]@series[["m3"]]@values
y_truth<-ts@series[["m3"]]@values
y_truth_reads<-ts_reads@series[["m3"]]@values

y_predict_rel <- exp(y_predict) / colSums(exp(y_predict))
y_truth_rel<-exp(y_truth)/colSums(exp(y_truth))
y_truth_rel[y_truth_reads==0] <- 0

y_predict_rel_3<-y_predict_rel
y_truth_rel_3<-y_truth_rel

##

fits[["hold-out-m2"]] <- mbtransfer(ts[c(2:4)], P = 2, Q = 2)
ts_preds[["hold-out-m2"]] <- predict(fits[["hold-out-m2"]], ts_missing[c(1)])

y_predict<-ts_preds[["hold-out-m2"]]@series[["m2"]]@values
y_truth<-ts@series[["m2"]]@values
y_truth_reads<-ts_reads@series[["m2"]]@values

y_predict_rel <- exp(y_predict) / colSums(exp(y_predict))
y_truth_rel<-exp(y_truth)/colSums(exp(y_truth))
y_truth_rel[y_truth_reads==0] <- 0

y_predict_rel_2<-y_predict_rel
y_truth_rel_2<-y_truth_rel

output_folder <- '../_outputs/log'

#dir.create(output_folder)

write.csv(y_truth_rel_5, paste(output_folder,"/m5_ground_truth.csv",sep=''))
write.csv(y_truth_rel_4, paste(output_folder,"/m4_ground_truth.csv",sep=''))
write.csv(y_truth_rel_3, paste(output_folder,"/m3_ground_truth.csv",sep=''))
write.csv(y_truth_rel_2, paste(output_folder,"/m2_ground_truth.csv",sep=''))

write.csv(y_predict_rel_5, paste(output_folder,"/m5_forecast.csv",sep=''))
write.csv(y_predict_rel_4, paste(output_folder,"/m4_forecast.csv",sep=''))
write.csv(y_predict_rel_3, paste(output_folder,"/m3_forecast.csv",sep=''))
write.csv(y_predict_rel_2, paste(output_folder,"/m2_forecast.csv",sep=''))


