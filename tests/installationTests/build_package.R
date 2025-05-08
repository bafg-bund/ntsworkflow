

# Build package ####

build_vignettes_to_inst <- function() {
  devtools::build_vignettes() # Builds vignettes to 'doc' and 'Meta'. Updates '.gitignore'.
  unlink(c("inst/doc", "inst/Meta"), recursive = TRUE) # Remove the directories if they exist
  dir.create("inst/doc"); dir.create("inst/Meta") # Create empty directories
  has_worked <- c( # Copy files to 'inst' subfolders
    file.copy(list.files("doc", full.names = TRUE), to = "inst/doc") 
    , file.copy(list.files("Meta", full.names = TRUE), to = "inst/Meta")
  )
  unlink(c("doc", "Meta"), recursive = TRUE) # Optional: Remove unwanted directories
  return(all(has_worked)) # Returns TRUE if everything worked OK
}

build_vignettes_to_inst() 

# only run this on linux
devtools::build(path = file.path(Sys.getenv("RENV_PATHS_CELLAR"), "ntsworkflow", glue::glue("ntsworkflow_{desc::desc_get_version()}.tar.gz")), binary = TRUE)

renv::deactivate()

install.packages("/beegfs/nts/renv_package_cellar/ntsworkflow/ntsworkflow_0.2.6.tar.gz")
<<<<<<< HEAD
library(ntsworkflow)
=======

>>>>>>> a131f2a (update to R 4.5.0)

renv::activate()
# only run this on windows
# delete .o and .so files in the src folder before building package in windows
# unlink(list.files("src", pattern = ".*\\.o$", f = TRUE))
# devtools::build(path = "Z:/G/G2/HRMS/sw_entwicklung/ntsworkflow/windows", binary = TRUE)
# unlink(c("src-i386", "src-x64"), recursive = TRUE)


#devtools::build(binary = TRUE)
#devtools::install(build_vignettes = TRUE)

# Copyright 2025 Bundesanstalt für Gewässerkunde
# This file is part of ntsworkflow

