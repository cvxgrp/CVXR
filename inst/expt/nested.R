setClass("L1", representation = list(x = "numeric", y = "character"))
setGeneric("setX<-", function(object, value) standardGeneric("setX<-") )
setMethod(f = "setX<-",
          signature = "L1",
          definition = function(object, value) {
              object@x <- value
              object
          })

setClass("L2", representation = list(l1 = "L1", z = "numeric"))
setGeneric("setL1<-", function(object, value) standardGeneric("setL1<-") )
setMethod(f = "setL1<-",
          signature = "L2",
          definition = function(object, value) {
              object@l1 <- value
              object
          })

##setClass("L3", representation = list(l2 = "L2"))

l1 <- new("L1", x=1:5, y=letters[1:5])
l2 <- new("L2", l1=l1, z=1:3)

l1p <- new("L1", x=10, y=letters[1:10])
l2p <- new("L2", l1=l1, z=4)

rm(l1, l1p)
for (bar

## Goal to replace X in l2$l1

