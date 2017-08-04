template <- '
<WHAT> = function(value) {
          if (missing(value)) {
              private$<WHAT> <- value
          } else {
              private$<WHAT>
          }
}
'

genActive <- function(what) {
    gsub("<WHAT>", what, template)
}

for (x in c("LP_CAPABLE", "SOCP_CAPABLE",
                            "SDP_CAPABLE",
                            "EXP_CAPABLE",
            "MIP_CAPABLE")) {
    catn(",")
    catn(genActive(x))
}

