library(RCurl)
url <- "ftp://ftp.bafg.de/pub/REFERATE/g2/ntsworkflow/"
userpwd <- "anonymous:"
filenames <- getURL(url, userpwd = userpwd,
                    ftp.use.epsv = FALSE, dirlistonly = TRUE) 
dat <- try(getURL("ftp://ftp.bafg.de/pub/REFERATE/g2/ntsworkflow/test2.csv", userpwd = userpwd))
dat
read.table(text = dat, sep = ",", header = T)
