% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility_functions.R
\name{scMikoUpdate}
\alias{scMikoUpdate}
\title{Reinstall scMiko package}
\usage{
scMikoUpdate(token = "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", ...)
}
\arguments{
\item{token}{autherization key to for private git repository}

\item{...}{Additional arguments passed to devtools::install_github()}
}
\description{
Reinstall scMiko package from private repository (https://github.com/NMikolajewicz/scMiko). Package is reloaded after update.
}
\examples{

# reinstall scMiko
scMikoUpdate()

}
\author{
Nicholas Mikolajewicz
}
