#' external_graphs
#'
#' Internal function to change rstudio default graph device
#' @param ext whether to switch on the external graph device
#' @keywords internal

#' @return a command
#' @examples \dontrun{
#' }
#'

external_graphs <- function(ext = TRUE) {

  # http://thecoatlessprofessor.com/programming/detecting-if-r-is-in-rstudio-
  # and-changing-rstudios-default-graphing-device/
  is.rstudio = function() {
    .Platform$GUI == "RStudio"
  }

  if( is.rstudio() ){
    if( isTRUE(ext) ){
      o = tolower(Sys.info()["sysname"])
      a = switch(o,
                 "darwin"  = "quartz",
                 "linux"   = "x11",
                 "windows" = "windows")
      options("device" = a)
    } else{
      options("device"="RStudioGD")
    }
    # Kill open graphic devices
    graphics.off()
  }
}


