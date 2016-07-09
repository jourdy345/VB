print.va <- function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
      if (length(x$tau)) {
        cat("\nQuantile:\n", x$tau,
        "\n\n")
    }
    if (length(coef(x))) {
        cat("Coefficients:\n")
        print.default(format(coef(x), digits = digits), print.gap = 2, 
            quote = FALSE)
    }else {cat("No linear effects\n")}
    cat("\n")
    cat("\nNumber of iterations:\n", x$i, 
        "\n\n", sep = "")
      if(length(which(names(x$gam)%in%c("geo", "random")==FALSE))){ 
         cat("Number of nonlinear effects:\n", length(which(names(x$gam)%in%c("geo", "random")==FALSE)),
         "\n\n")
     } 
       if(length(which(names(x$gam)%in%c("geo")))){ 
         cat("Spatial effect for :", length(x$gam[[which(names(x$gam)%in%c("geo"))]]), "locations",
         "\n\n")
     }
       if(length(which(names(x$gam)%in%c("random")))){ 
         cat("Random effects for", length(x$gam[[which(names(x$gam)%in%c("random"))]]),"individuals/ groups",
         "\n\n")
     }#else cat("No nonlinear effects\n")
#     cat("\n")
    invisible(x)
   #  NextMethod("print")
}
