safeguarding <- function(.f) {
  tryCatch(
    .f,
    warning = function(w)
      NULL,
    error = function(e)
      NULL
  )
}
