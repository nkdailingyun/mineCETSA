.onLoad <- function(...) {

  message <- paste("Welcome to mineCETSA!",
    "Now it's the time to start mining your CETSA data!",
    "The two vignettes below explain in detail how to use this package",
    "Access the respective vignettes by key in: ",
    "vignette('CETSA_melt_curve', package='mineCETSA')", "or",
    "vignette('CETSA_ITDR_ITTR', package='mineCETSA')", "or",
    "vignette('CETSA_2D', package='mineCETSA')", sep="\n")

  packageStartupMessage(message)

}
