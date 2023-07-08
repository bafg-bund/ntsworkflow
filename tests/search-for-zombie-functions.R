
library(glue)

tosearch <- "FindPeaks_BfG_scripting"

system(glue("grep -rn '{tosearch}' tests"))
system(glue("grep -rn '{tosearch}' R"))
system(glue("grep -rn '{tosearch}' inst/shiny"))
system(glue("grep -rn '{tosearch}' ~/projects/ntsautoeval"))
