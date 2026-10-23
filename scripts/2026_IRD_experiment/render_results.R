# Refit and regenerate the consolidated results, prediction tables and PDF.
# All analysis code is kept in outputs/2026_IRD_experiment/wheat_competition_results.Rmd.
arg <- grep('^--file=',commandArgs(),value=TRUE)
root <- dirname(dirname(dirname(normalizePath(sub('^--file=','',arg[1])))))
pandoc <- 'C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools'
if (!rmarkdown::pandoc_available() && file.exists(file.path(pandoc,'pandoc.exe')))
  Sys.setenv(RSTUDIO_PANDOC=pandoc)
rmarkdown::render(file.path(root,'outputs','2026_IRD_experiment',
  'wheat_competition_results.Rmd'),clean=TRUE)