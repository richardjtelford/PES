fortify.envfit <- function (x, choices = c(1, 2), labels, at = c(0, 0), axis = FALSE, p.max = NULL, ...) 
{
  formals(arrows) <- c(formals(arrows), alist(... = ))
  labs <- list(v = rownames(x$vectors$arrows), f = rownames(x$factors$centroids))
  if (!missing(labels)) {
    if (is.list(labels)) {
      if (!is.null(labs$v) && !is.null(labels$vectors)) 
        labs$v <- labels$vectors
      if (!is.null(labs$f) && !is.null(labels$factors)) 
        labs$f <- labels$factors
    }
    else {
      if (!is.null(x$vectors) && !is.null(x$factors)) 
        stop("needs a list with both 'vectors' and 'factors' labels")
      if (!is.null(x$factors)) 
        labs$f <- labels
      else labs$v <- labels
    }
  }
  vect <- NULL
  if (!is.null(p.max)) {
    if (!is.null(x$vectors)) {
      take <- x$vectors$pvals <= p.max
      x$vectors$arrows <- x$vectors$arrows[take, , drop = FALSE]
      labs$v <- labs$v[take]
      x$vectors$r <- x$vectors$r[take]
      if (nrow(x$vectors$arrows) == 0) 
        x$vectors <- vect <- NULL
    }
    if (!is.null(x$factors)) {
      tmp <- x$factors$pvals <= p.max
      nam <- names(tmp)[tmp]
      take <- x$factors$var.id %in% nam
      x$factors$centroids <- x$factors$centroids[take, , drop = FALSE]
      labs$f <- labs$f[take]
      if (nrow(x$factors$centroids) == 0) 
        x$factors <- NULL
    }
  }
  if (!is.null(x$vectors)) {
    vect <- sqrt(x$vectors$r) * x$vectors$arrows[, choices, drop = FALSE]
    if (axis) {
      maxarr <- round(sqrt(max(x$vectors$r)), 1)
      ax <- -c(-1, 0, 1) * maxarr
    }
    vect <- arrow.mul * vect
    vect <- sweep(vect, 2, at, "+")
  }
  data.frame(vect, labs = labs$v)
}
