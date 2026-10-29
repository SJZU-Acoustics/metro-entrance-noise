suppressPackageStartupMessages(library(readr))
trace("read_csv", where=asNamespace("readr"), print=FALSE, tracer=quote({
 if(is.character(file) && length(file)==1 && file.exists(file)) cat(normalizePath(file), "\n",file="verification/input_reads.txt",append=TRUE)
}))
