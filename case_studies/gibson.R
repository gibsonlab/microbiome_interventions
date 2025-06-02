#install.packages("devtools")
#install.packages("tidyverse")
#devtools::install_github("gathanei/xyz")
#devtools::install_github("krisrs1128/mbtransfer")


library(mbtransfer)
library(tidyverse)
library(glue)


subject <- read_csv("../_data/healthy/subject.csv")
interventions <- read_csv("../_data/healthy/interventions.csv") |> column_to_rownames("sample")
reads <- read_csv("../_data/healthy/reads.csv") |> column_to_rownames("sample")
samples <- read_csv("../_data/healthy/samples.csv")

eps<-1e-4

R<-as.matrix(reads)
R<- R / rowSums(R)
R <- log10(R+eps)

ts <- R |>
  ts_from_dfs(interventions, samples, subject) |>
  interpolate(method = "linear")

ts

fits <- list()
ts_preds <- list()
subject_data(ts)

# mice are m2, m3, m4, m5 ... so  here we m2,m3,m4 and then tests on m5, starting from only knowing the initial condition

ts_missing <- subset_values(ts, 1)

fits[["hold-out-m5"]] <- mbtransfer(ts[c(1:3)], P = 2, Q = 2)
ts_preds[["hold-out-m5"]] <- predict(fits[["hold-out-m5"]], ts_missing[c(4)])

fits[["hold-out-m4"]] <- mbtransfer(ts[c(1:2, 4)], P = 2, Q = 2)
ts_preds[["hold-out-m4"]] <- predict(fits[["hold-out-m4"]], ts_missing[c(3)])

fits[["hold-out-m3"]] <- mbtransfer(ts[c(1, 3:4)], P = 2, Q = 2)
ts_preds[["hold-out-m3"]] <- predict(fits[["hold-out-m3"]], ts_missing[c(2)])

fits[["hold-out-m2"]] <- mbtransfer(ts[c(2:4)], P = 2, Q = 2)
ts_preds[["hold-out-m2"]] <- predict(fits[["hold-out-m2"]], ts_missing[c(1)])


df <- map_dfr(ts_preds, ~ reshape_preds(ts, .), .id = "generalization")
rmse <- sqrt(mean((df$y_hat - df$y)^2, na.rm = TRUE))
print(rmse)

# Calculate RMSE for subject "m5"
rmse <- function(actual, predicted) {
  sqrt(mean((actual - predicted)^2, na.rm = TRUE))
}

rmse_m5 <- df %>%
  filter(subject == "m5") %>%
  summarize(rmse_value = rmse(y_hat, y))

print(rmse_m5)

y1<-ts_preds[["hold-out-m5"]]@series[["m5"]]@values
y2<-ts@series[["m5"]]@values
out<-rmse(y1, y2)
print(out)

# only calculate error for the taxa-timepoints that have nonzero reads like we do in the analysis in our paper
# technically we should not compare to "ground truth" at timepoints that are interpolated either but havent done this yet
# we should remove the first time point as well

y3<-y2
y3[y3 == log10(eps)] <- NA #mask values that originally had zero reads 
out<-rmse(y1, y3)
print(out)


output_folder <- '../_outputs'

write.csv(ts@series[["m5"]]@values, paste(output_folder,"/m5_ground_truth.csv",sep=''))
write.csv(ts@series[["m4"]]@values, paste(output_folder,"/m4_ground_truth.csv",sep=''))
write.csv(ts@series[["m3"]]@values, paste(output_folder,"/m3_ground_truth.csv",sep=''))
write.csv(ts@series[["m2"]]@values, paste(output_folder,"/m2_ground_truth.csv",sep=''))

write.csv(ts_preds[["hold-out-m5"]]@series[["m5"]]@values, paste(output_folder,"/m5_forecast.csv",sep=''))
write.csv(ts_preds[["hold-out-m4"]]@series[["m4"]]@values, paste(output_folder,"/m4_forecast.csv",sep=''))
write.csv(ts_preds[["hold-out-m3"]]@series[["m3"]]@values, paste(output_folder,"/m3_forecast.csv",sep=''))
write.csv(ts_preds[["hold-out-m2"]]@series[["m2"]]@values, paste(output_folder,"/m2_forecast.csv",sep=''))


