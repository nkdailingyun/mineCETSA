.onLoad <- function(...) {

  message <- paste("The vignettes below explain in detail how to use mineCETSA package",
               "Access the respectie vignettes by key in:",
               "vignette('CETSA_melt_curve', package='mineCETSA')", "or",
               "vignette('CETSA_ITDR_ITTR', package='mineCETSA')", sep="\n")

  packageStartupMessage(message)

}
