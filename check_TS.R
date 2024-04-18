
## function to check multivariate TS object to make sure it is in the format
## we expect for Lund workshop
 check_TS <- function (mvTS){
   
   msg <- character() # vector to hold warning messages
   
   # check if nsamp dimensions are correct
   if(mvTS$nsamp != nrow(mvTS$M)) msg <- append(msg, "# of rows of M does not match nsamp.")
   if(mvTS$nsamp != length(mvTS$S)) msg <- append(msg, "length of S does not match nsamp.")
   if(mvTS$nsamp != length(mvTS$nn)) msg <- append(msg, "length of nn does not match nsamp.")
   if(mvTS$nsamp != length(mvTS$tt)) msg <- append(msg, "length of tt does not match nsamp.")
   
   # check if nvar dimensions are correct
   if(mvTS$nvar != ncol(mvTS$M)) msg <- append(msg, "# of columns of M does not match nvar.")
   Sdim <- sapply(mvTS$S, dim)
   if(!all(Sdim == mvTS$nvar)) msg <- append(msg, "dimensions of S matrices are not all nvar x nvar.")
   
   # check if tt starts at zero and strictly increases; also check units
   if(mvTS$tt[1] != 0) msg <- append(msg, "tt of the first sample is not zero.")
   dtt <- diff(mvTS$tt)
   if(any(dtt <= 0)) msg <- append(msg, "tt is not strictly increasing.")
   if(!mvTS$time.units %in% c("Myr", "yr")) msg <- append(msg, "time.units is not either 'Myr' or 'yr'.")
   
   # check if there are any NA's anywhere
   checkNA <- sapply(mvTS, anyNA)
   if(any(checkNA)) {
     NAmsg <- paste0("These elements had NA's: ", names(checkNA)[checkNA])
     msg <- append(msg, NAmsg) 
   }

   # check if other informational elements are included
   if(is.null(mvTS$taxon)) msg <- append(msg, "'taxon' not included.") 
   if(is.null(mvTS$sex)) msg <- append(msg, "'sex' not included.")
   if(is.null(mvTS$reference)) msg <- append(msg, "'reference' not included.")
   #if(is.null(mvTS$transformation)) msg <- append(msg, "'transoformation' not included.")
   #if(is.null(mvTS$trait.dim)) msg <- append(msg, "'trait.dim' not included.")
   
   
   if(length(msg) == 0) cat("The mvTS object looks good.\n") else{
     cat("Note the following warnings:", "\n")
     for(i in 1:length(msg))
       cat("  ", msg[i], "\n")
     cat("\n")
    }
   
   invisible(msg)
 }
 
 
