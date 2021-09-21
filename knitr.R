setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

knit(input="readme.rmd", output = "readme.md")
knit(input="LAB1/readme.rmd", output = "LAB1/readme.md")
knit(input="LAB2/readme.rmd", output = "LAB2/readme.md")
