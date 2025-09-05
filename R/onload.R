#' Package startup message
#'
#' This function is called when the package is loaded and displays
#' information about the package creators.
#'
#' @param libname Library path
#' @param pkgname Package name
#' @export
.onLoad <- function(libname, pkgname) {
  # Get package description
  desc <- packageDescription(pkgname, lib.loc = libname)

  # Extract version and author information
  version <- desc$Version
  authors <- paste(desc$Author, collapse = " and ")

  # Create and display the message
  message <- paste0("manhattantwin ", version, " created by ", authors)
  packageStartupMessage(message)
}
