##
## We shadow a CPP class with an R6 class using the (almost) same class name!
## Some utility functions
##
rcppMungedName <- function(cppClassName, methodName, thisPkg = getPackageName()) {
    ## Remove the CVXCanon. prefix from all the class names.
    paste0(thisPkg, "_", gsub("^CVXCanon\\.", "", cppClassName) , "__", methodName)
}

