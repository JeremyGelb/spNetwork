% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{esc_kernel_loo}
\alias{esc_kernel_loo}
\title{The recursive function to calculate continuous NKDE likelihood cv}
\arguments{
\item{kernel_func}{a cpp pointer function (selected with the kernel name)}

\item{samples_k}{a numeric vector of the actual kernel values, updated at
each recursion}

\item{neighbour_list}{a List, giving for each node an IntegerVector with
its neighbours}

\item{edge_mat}{matrix, to find the id of each edge given two neighbours.}

\item{v}{the actual node to consider for the recursion (int)}

\item{v1}{the connected node to consider for the recursion (int)}

\item{l1}{the edge connecting v and v1 (int)}

\item{d}{the actual distance traveled before the recursion}

\item{alpha}{the actual alpha value before the recursion}

\item{bw}{the kernel bandwidth}

\item{line_weights}{a vector with the length of the edges}

\item{depth}{the actual recursion depth}

\item{max_depth}{the maximum recursion depth}
}
\value{
a vector with the kernel values calculated for each samples from
the first node given
}
\description{
The recursive function to calculate continuous NKDE likelihood cv
}
