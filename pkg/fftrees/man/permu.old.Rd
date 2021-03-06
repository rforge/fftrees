% Generated by roxygen2 (4.0.2): do not edit by hand
\name{permu.old}
\alias{permu.old}
\title{permu.old}
\usage{
permu.old(perm, fun, steps = FALSE, current = NULL)
}
\arguments{
\item{perm}{numeric vector; numbers which should be permuted}

\item{fun}{function; callback function with 1 parameter}

\item{steps}{logical; Generate only full permutations, or generate steps between them? Default: \code{FALSE}}

\item{current}{numeric vector; Prepends given vector to permutation before calling \code{fun}}
}
\value{
nothing. If you want to save calculations, you need to write them into global workspace via \code{fun}
}
\description{
Generates permutations of \code{perm} and calls a callback \code{fun} for every permutation, without preallocating
memory. This helps, when permutating very large vectors.
}

