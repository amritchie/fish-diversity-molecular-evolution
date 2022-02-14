my.read.GenBank <- function (access.nb, seq.names = access.nb, species.names = TRUE, 
          gene.names = FALSE, as.character = FALSE) 
{
  N <- length(access.nb)
  a <- 1L
  b <- if (N > 200) 
    200L
  else N
  fl <- tempfile()
  repeat {
    URL <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=", 
                  paste(access.nb[a:b], collapse = ","), "&rettype=fasta&retmode=text")
    
    X <- scan(file = URL, what = "", sep = "\n", quiet = TRUE)
    cat(X, sep = "\n", file = fl, append = TRUE)
    if (b == N) 
      break
    a <- b + 1L
    b <- b + 200L
    if (b > N) 
      b <- N
  }
  res <- read.FASTA(fl)
  if (is.null(res)) 
    return(NULL)
  attr(res, "description") <- names(res)
  if (length(access.nb) != length(res)) {
    names(res) <- gsub("\\..*$", "", names(res))
    failed <- paste(access.nb[!access.nb %in% names(res)], 
                    collapse = ", ")
    warning(paste0("cannot get the following sequences:\n", 
                   failed))
  }
  else names(res) <- access.nb
  if (as.character) 
    res <- as.character(res)
  if (species.names) {
    a <- 1L
    b <- if (N > 200) 
      200L
    else N
    sp <- character(0)
    repeat {
      URL <- paste("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=", 
                   paste(access.nb[a:b], collapse = ","), 
                   "&rettype=gb&retmode=text", sep = "")
      X <- scan(file = URL, what = "", sep = "\n", 
                quiet = TRUE, n = -1)
      sp <- c(sp, gsub(" +ORGANISM +", "", 
                       grep("ORGANISM", X, value = TRUE)))
      if (b == N) 
        break
      a <- b + 1L
      b <- b + 200L
      if (b > N) 
        b <- N
    }
    attr(res, "species") <- gsub(" ", "_", 
                                 sp)
  }
  if (gene.names) 
    warning("you used 'gene.names = TRUE': this option is obsolete; please update your code.")
  res
}