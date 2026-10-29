suppressPackageStartupMessages(library(readr))
x<-read_csv("intermediate/collocation_pairs.csv",show_col_types=FALSE)
d<-x$recorder_db-x$reference_db
stopifnot(nrow(x)==14,all(is.finite(d)))
write_csv(data.frame(n=length(d),mean_bias_db=mean(d),rmse_db=sqrt(mean(d^2))),"intermediate/analysis/collocation_summary.csv")
