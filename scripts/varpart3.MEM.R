varpart3.MEM <- 
	function(Y, X1, X2, X3, is.MEM=NULL, method="hierarchical")
#
# Hierarchical and proportional variation partitioning of Y 
# with respect to X1, X2 and X3.
#
# Usage:
# Y : Vector, matrix or data frame of response data.
# X1, X2, X3 : Vectors, matrices or data frames of explanatory data.
# is.MEM=NULL : No matrix of MEM eigenfunctions, or only one. That matrix 
#       will receive no special treatment since no hierarchy of spatial 
#       submodels has to be set. Classical variation partitioning is used.
# is.MEM = c(2,3) : Matrices X2 and X3 contain MEM eigenfunctions. Place the 
#		broad-scale MEM model before the fine-scale model.
# is.MEM = c(1:3) : Matrices X1, X2 and X3 contain MEM eigenfunctions. Place the 
#		broad-scale MEM model before the midle-scale model, which must be placed 
#		before the fine-scale model.
# method="hierarchical" : Hierarchical partitioning.
#       ="proportional" : Proportional apportioning of the shared fractions.
#       ="classical"    : classical variation partitioning.
# Methods can be abbreviated to any recognizable subset, e.g. "h", "p", "c".
#
# Note -- Classical variation partitioning is used if method="classical",
# if is.MEM=NULL, or if there is a single file of MEM (e.g. is.MEM=3).
#
# Author:: Pierre Legendre
# License: GPL-2
{
require(vegan)
method <- match.arg(method, c("hierarchical", "proportional","classical"))
x <- varpart(Y, X1, X2, X3)   # The heavy work is done by varpart()
#
	if(method=="classical") return(x)
	if(any(is.MEM < 1) || any(is.MEM > 3)) stop("Error in is.MEM statement")
    if(!is.null(is.MEM)) {
    	if(length(is.MEM) == 1) {
    	   cat("Partitioning with no special treatment for a single MEM file\n")
    	   return(x)
    	   }
    	if(length(is.MEM) > 1) {
    		are.MEM <- TRUE
    		# cat("length(is.MEM) =", length(is.MEM),"\n")
    		# cat("(length(is.MEM) == 2) =", (length(is.MEM) == 2),"\n")
    		if(length(is.MEM) == 2) {
    			# cat("sum(is.MEM - c(2,3)) =", sum(is.MEM - c(2,3)),"\n")
    			if(sum(is.MEM - c(2,3)) != 0) {
    				stop("The MEM files must be placed last")
    				} else {
    				cat("Files X2 and X3 contain MEM eigenfunctions\n") }
    			}
    		if(length(is.MEM) == 3) {
    			cat("Files X1, X2 and X3 contain MEM eigenfunctions\n") }
    		}
    	} else { return(x) }
#
if(length(is.MEM) == 2) {    # Two matrices of MEM
	### method="hierarchical"
	if(method=="hierarchical") {
		cat("Hierarchical partitioning of the shared fractions\n")
		x$part$indfract[2,3] <- sum(x$part$indfract[c(2,5),3])
		x$part$indfract[4,3] <- sum(x$part$indfract[c(4,7),3])
		x$part$indfract[c(5,7),3] <- NA
		} else {
	### method="proportional"
		cat("Proportional apportioning of the shared fractions\n")
		bd <- x$part$contr1[3,3]
		cf <- x$part$contr1[6,3]
		b <- x$part$indfract[2,3]
		c <- x$part$indfract[3,3]
		d <- x$part$indfract[4,3]
		e <- x$part$indfract[5,3]
		f <- x$part$indfract[6,3]
		g <- x$part$indfract[7,3]
		e.to.b <- e * bd/(bd+cf)
		e.to.c <- e * cf/(bd+cf)
		g.to.d <- g * bd/(bd+cf)
		g.to.f <- g * cf/(bd+cf)
		#
		x$part$indfract[2,3] <- b + e.to.b
		x$part$indfract[3,3] <- c + e.to.c
		x$part$indfract[4,3] <- d + g.to.d
		x$part$indfract[6,3] <- f + g.to.f
		x$part$indfract[c(5,7),3] <- NA
		}
	# Check results
	if(abs(sum(x$part$indfract[c(2:4,6),3])-x$part$fract[6,3]) > 1e-10) 
		cat("Error in fractions for X2+X3\n")
	return(x)
	
	} else if(length(is.MEM) == 3) {  # Three matrices of MEM
	### method="hierarchical"
	if(method=="hierarchical") {
		cat("Hierarchical partitioning of the shared fractions\n")
		x$part$indfract[1,3] <- sum(x$part$indfract[c(1,4,6,7),3])
		x$part$indfract[2,3] <- sum(x$part$indfract[c(2,5),3])
		x$part$indfract[4:7,3] <- NA
		} else {
	### method="proportional"
		cat("Proportional apportioning of the shared fractions\n")
		a <- x$part$indfract[1,3]
		b <- x$part$indfract[2,3]
		c <- x$part$indfract[3,3]
		d <- x$part$indfract[4,3]
		e <- x$part$indfract[5,3]
		f <- x$part$indfract[6,3]
		g <- x$part$indfract[7,3]
		d.to.a <- d * a/(a+b)
		d.to.b <- d * b/(a+b)
		e.to.b <- e * b/(b+c)
		e.to.c <- e * c/(b+c)
		f.to.a <- f * a/(a+c)
		f.to.c <- f * c/(a+c)
		g.to.a <- g * a/(a+b+c)
		g.to.b <- g * b/(a+b+c)
		g.to.c <- g * c/(a+b+c)
		#
		x$part$indfract[1,3] <- a + d.to.a + f.to.a + g.to.a
		x$part$indfract[2,3] <- b + d.to.b + e.to.b + g.to.b
		x$part$indfract[3,3] <- c + e.to.c + f.to.c + g.to.c
		x$part$indfract[4:7,3] <- NA
		}
	# Check results
	if(abs(sum(x$part$indfract[1:3,3])-x$part$fract[7,3]) > 1e-10) 
		cat("Error in fractions for X1+X2+X3\n")
	return(x)
	}
}