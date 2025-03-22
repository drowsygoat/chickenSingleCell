#' List Objects with Sizes (Shortened Name: los)
#'
#' This function `los` lists all objects in the specified environment, sorted by their memory sizes in a human-readable format. It allows for quick identification of large objects that may be using significant memory.
#'
#' @param envir Environment to search for objects, defaults to global environment.
#' @param order Descending or ascending order for size sorting. Default is 'decreasing'.
#' @return Prints a table of object names and their formatted sizes. Does not return any value.
#' @examples
#' los()
#' los(envir = .GlobalEnv, order = "decreasing")
#' @export
los <- function(envir = .GlobalEnv, order = "decreasing") {
  # Input validation using stopifnot
  stopifnot(inherits(envir, "environment"), order %in% c("increasing", "decreasing"))

  # Get the list of objects in the specified environment
  objects <- ls(envir = envir)
  # Calculate sizes and format them
  sizes <- sapply(objects, function(x) object.size(get(x, envir = envir)), simplify = FALSE)
  formatted_sizes <- sapply(sizes, function(x) format(x, units = "auto"))
  # Sort the list of sizes based on the specified order
  object_sizes <- order(sapply(objects, function(x) object.size(get(x, envir = envir))), decreasing = order == "decreasing")
  
  # Create a data frame for pretty printing
  size_data <- data.frame(
    Object = objects[object_sizes],
    Size = formatted_sizes[object_sizes],
    stringsAsFactors = FALSE
  )

  # Print the data frame
  print(size_data, row.names = FALSE)
}