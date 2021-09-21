setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

knit(input="readme.rmd", output = "readme.md")
