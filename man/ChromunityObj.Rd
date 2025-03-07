% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chromunity.R
\name{ChromunityObj}
\alias{ChromunityObj}
\title{Chromunity object}
\description{
Vectorized object that stores the output of chromunity analysis and can be inputted into annotate / synergy functions.
The main accessors are $binsets and $concatemers which each return GRanges linked through a chromunity id $chid field

Chromunities can be subsetted, concatenated.  Concatenation will result in deduping any $chid that are present in two inputs.
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-Chromunity-new}{\code{ChromunityObj$new()}}
\item \href{#method-Chromunity-subset}{\code{ChromunityObj$subset()}}
\item \href{#method-Chromunity-print}{\code{ChromunityObj$print()}}
\item \href{#method-Chromunity-gtrack}{\code{ChromunityObj$gtrack()}}
\item \href{#method-Chromunity-gm}{\code{ChromunityObj$gm()}}
\item \href{#method-Chromunity-clone}{\code{ChromunityObj$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Chromunity-new"></a>}}
\if{latex}{\out{\hypertarget{method-Chromunity-new}{}}}
\subsection{Method \code{new()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ChromunityObj$new(
  binsets = GRanges(),
  concatemers = GRanges(),
  meta = data.table(),
  verbose = TRUE
)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Chromunity-subset"></a>}}
\if{latex}{\out{\hypertarget{method-Chromunity-subset}{}}}
\subsection{Method \code{subset()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ChromunityObj$subset(ix, ...)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Chromunity-print"></a>}}
\if{latex}{\out{\hypertarget{method-Chromunity-print}{}}}
\subsection{Method \code{print()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ChromunityObj$print(...)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Chromunity-gtrack"></a>}}
\if{latex}{\out{\hypertarget{method-Chromunity-gtrack}{}}}
\subsection{Method \code{gtrack()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ChromunityObj$gtrack(
  name = "",
  binsets = TRUE,
  concatemers = TRUE,
  heatmap = FALSE,
  binsize = 10000,
  clim = NA,
  ...
)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Chromunity-gm"></a>}}
\if{latex}{\out{\hypertarget{method-Chromunity-gm}{}}}
\subsection{Method \code{gm()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ChromunityObj$gm(binsize = 5000)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-Chromunity-clone"></a>}}
\if{latex}{\out{\hypertarget{method-Chromunity-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{ChromunityObj$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
