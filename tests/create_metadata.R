
# Script to create a metadata file
# Will have the same name and path as RDS file
# usage Rscript create_metadata.R path_to_RDS -r Rhein -lon 7.43232 -lat 52.09283 -s rhein_ko_l -m Wasser -g 2
ca <- commandArgs(TRUE)
filePath <- ca[1]
riv <- "river_name"
if (any(grepl("-r", ca)))
  riv <- ca[grep("-r", ca) + 1]

lon <- 0.1
if (any(grepl("-lon", ca)))
  lon <- as.numeric(ca[grep("-lon", ca) + 1])
stopifnot(!is.na(lon))

lat <- 0.1
if (any(grepl("-lat", ca)))
  lat <- as.numeric(ca[grep("-lat", ca) + 1])
stopifnot(!is.na(lat))

station <- "station_name"
if (any(grepl("-s", ca)))
  station <- ca[grep("-s", ca) + 1]

matrix <- "matrix_name"
if (any(grepl("-m", ca)))
  matrix <- ca[grep("-m", ca) + 1]

gtk <- 0
if (any(grepl("-g", ca)))
  gtk <- as.numeric(ca[grep("-g", ca) + 1])
stopifnot(!is.na(gtk))

create_meta_csv <- function(path, filename, start, lon, lat, station, matrix, river, gkz) {
  df <- data.frame(filename, start, station, lon, lat, matrix, river, gkz)
  write.csv(df, file = path, row.names = F)
}
filePath <- normalizePath(filePath)

stopifnot(file.exists(filePath))
if (!grepl("\\.RDS$", filePath))
  stop("That is not an RDS file")
all_data <- readRDS(filePath)
samp <- subset(all_data$sampleList, !deleted)

stopifnot(is.data.frame(samp))
filePath2 <- sub("\\.RDS$", ".csv", filePath)

create_meta_csv(
  filePath2, 
  basename(samp$File),
  lubridate::ymd(stringr::str_extract(basename(samp$File), "\\d{8}")),
  lon,
  lat,
  station,
  matrix,
  riv,
  gtk
)

message("done")
