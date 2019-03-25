\name{baggedImpute}
\alias{baggedImpute}

\title{
Implement bootstrap aggregation for the imputation of missing values in time series
}
\description{
A function that applies non-overlapping block bootstrap (NBB) (Carlstein, 1986), overlapping block bootstrap (MBB) (Kunsch, 1989; Liu and Singh, 1992), and cicular block bootstrap (CBB) (Politis and Romano, 1991) along with marked point processes (Parzen, 1963; Niebuhr, 2017) and bootstrap aggregating technique (Breiman, 1996) to impute missing values in time series.
}
\usage{
baggedImpute(data, bootrep = 1000, blocklen = 2, 
        BBM = c("MBB", "NBB", "CBB"), misstr = c("ams", "ivs"), 
        interpolation = c("Stineman", "Kalman", "WMA") )         
}
\arguments{
  \item{data}{
    An \code{array} including the time series with missing values. 
  }
  \item{bootrep}{
    The number of bootstrap replications. 
  }
  \item{blocklen}{
    Block length of NBB, MBB, or CBB.
  }
  \item{BBM}{
    Specifies the block bootstrap method. Takes one of \code{"MBB"}, \code{"NBB"}, or \code{"CBB"}.
  }
  \item{misstr}{
    Specifies the marked point process. Takes \code{"ams"} for amplitude modulated sequence and \code{"ivs"} integer-valued sequence. 
  }  
  \item{interpolation}{
   Specifies the the interpolation method to be applied.
  }
}
\details{

The details of the implementaion can be found in Andiojaya and Demirhan (2019).
}

\value{

\item{RegularSeries}{An object including the series with imputed missing values.}

}

\references{

Andiojaya, A., Demirhan, H. (2019). A bagging algorithm for the imputation of missing values in time series. Submitted to Expert Systems with Applications.

Breiman, L. (1996). Bagging predictors. Machine learning, 24(2):123-140.

Carlstein, E. (1986). The use of subseries values for estimating the variance of a general statistic from a stationary sequence. The Annals of Statistics, pages 1171-1179.

Kunsch, H. R. (1989). The jackknife and the bootstrap for general stationary observations. The Annals of Statistics, pages 1217-1241.

Liu,  R.  Y.  and  Singh,  K.  (1992).   Moving  blocks  jackknife  and  bootstrap capture weak dependence. InIn:  R. Lepage and L. Billard, Eds., Exploring the limits of bootstrap, 225-248. John Wiley, New York.

Niebuhr, T. (2017). Irregularly observed time series - some asymptotics and the block bootstrap. Statistics, 51(5):1118-1131.

Parzen, E. (1963). On spectral analysis with missing observations and amplitude modulation. Sankhya: The Indian Journal of Statistics, Series A, pages 38-392.

Politis, D. N. and Romano, J. P. (1991). A circular block-resampling procedure for stationary data. Technical Report 370, Stanford University, Department of Statistics.

}

\author{
Haydar Demirhan and Agung Andiojaya

Maintainer: Haydar Demirhan <haydar.demirhan@rmit.edu.au>
}

\examples{

data <- ts(rnorm(100, 0, 1))
data[ceiling(runif(10, 1,100))] <- NA

baggedImpute(data = data, bootrep = 10, blocklen = 2, 
             BBM = "NBB", misstr = "ams", 
             interpolation = "Stineman" ) 

}
