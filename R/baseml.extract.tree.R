require (fs)

baseml.extract.tree <- function (path)
{

    paml <- readLines(path)
    rex <- regmatches(paml, regexec("^(.*);$", paml))
    tr <- read.tree(text = rev(rex[sapply(rex, function (x) length(x) > 0)])[[1]][1])

    return (tr)

}

codeml.extract.trees <- function (path)
{

    paml <- readLines(path)
    rex <- regmatches(paml, regexec("^(.*);$", paml))
    trs <- sapply(1:3, function (x) read.tree(text = rev(
                                                     rex[sapply(rex, function (x) length(x) > 0)]
                                                 )[[x]][1]), simplify = F)
    names(trs) <- c("w", "dN", "dS")
    return (trs)

}
