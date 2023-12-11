# generates setGeneric saving typing
gen <- function(name, def) sprintf('setGeneric("%s", %s { standardGeneric("%s") })', name, def, name)
