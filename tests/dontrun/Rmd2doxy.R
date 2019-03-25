# conversion of Rmd formatting to Doxygen formatting

Rmd_lines <- readLines("markdown_math.tex")
# inline math
doxy_lines <- gsub("\\w*(?<!\\\\)[$]", "\\\\f$",
                   Rmd_lines, perl = TRUE)
# equation math
doxy_lines <- gsub("\\\\\\[", "\\\\f[", doxy_lines, perl = TRUE)
doxy_lines <- gsub("\\\\\\]", "\\\\f]", doxy_lines, perl = TRUE)
# unlink certain words
doxy_lines <- gsub("(Wishart)", "%\\1", doxy_lines, perl = TRUE)
# prepend doxygen comment
doxy_lines <- gsub("^", "/// ", doxy_lines, perl = TRUE)
# display output
cat(doxy_lines, sep = "\n")
