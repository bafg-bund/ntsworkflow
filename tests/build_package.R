# Build package


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

build_vignettes_to_inst() # Call the function

# only run this on linux
devtools::build(path = "~/G2_Z/HRMS/sw_entwicklung/ntsworkflow/linux", binary = TRUE)
# to copy the whole library to this directory as a tarball, use the following:
# first clean and rebuild!
# system("rm -r ~/R/myR36lib/00LOCK-ntsworkflow")
# system("tar czf ~/HRMS_Z/sw_entwicklung/ntsworkflow/linux/R36lib_ntsworkflow_0.2.0.tar.gz -C ~/R/myR36lib .")

# only run this on windows
# delete .o and .so files in the src folder before building package in windows
# unlink(list.files("src", pattern = ".*\\.o$", f = TRUE))
# devtools::build(path = "Z:/G/G2/HRMS/sw_entwicklung/ntsworkflow/windows", binary = TRUE)
# unlink(c("src-i386", "src-x64"), recursive = TRUE)


#devtools::build(binary = TRUE)
#devtools::install(build_vignettes = TRUE)

#just a test
