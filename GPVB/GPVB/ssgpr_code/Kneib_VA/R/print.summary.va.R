print.summary.va <- function(x, digits = max(3, getOption("digits") - 3), ...){
     cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
      if (length(x$tau)) {
        cat("\nQuantile:\n", x$tau,
        "\n\n")
    }
     if(length(x$coef)){
        cat("Coefficients:\n")
        print.default(format(x$coef, digits = digits), print.gap = 2, 
            quote = FALSE)
        cat( "\n\n")
       
     }
     if(length(x$nl)){
       cat("Nonlinear effects:\n",
       "there are",length(x$nl),"nonlinear effects:", paste(names(x$nl)), "\n\n")
     }
 
     if(length(x$geo)){
       cat("Spatial effect:\n",
       "the spatial variable contains ",length(x$geo[[2]]),
       "regions, the estimation was based on ",attr(x$geo, "number"),"regions", "\n\n")
     }
          if(length(x$random)){
       cat("Random effect:\n",
       "the random effect contains",length(x$random), "groups/individuals", "\n\n")
     }
 
}