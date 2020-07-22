safeguarding <- function(.f) {
  tryCatch(
    .f,
    error = function(e)
      NULL
  )
}
