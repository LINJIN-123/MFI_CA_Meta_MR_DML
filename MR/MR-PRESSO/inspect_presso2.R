library(MRPRESSO)
b <- deparse(body(mr_presso))
cat("total lines:", length(b), "\n\n")
cat(b[60:length(b)], sep = "\n")
