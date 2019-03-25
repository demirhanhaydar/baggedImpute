#' @import stats
#' @importFrom imputeTS na.interpolation
#' @importFrom imputeTS na.kalman
#' @importFrom imputeTS na.ma
#' @export

baggedImpute <- function(data, bootrep = 1000, blocklen = 2, BBM = c("MBB", "NBB", "CBB"), misstr = c("ams", "ivs"), interpolation = c("Stineman", "Kalman", "WMA") ){
  if (BBM == "CBB"){
    BootFunc <- c(GenCBB)  
  } else if (BBM == "MBB"){
    BootFunc <- c(GenMBB)  
  } else if (BBM == "NBB"){
    BootFunc <- c(GenNBB)  
  }
  listNA <- which(is.na(data))
  imputedSeries <- data
  if(misstr == 'ams'){    #combine with ams
    
    boot.result <- BootFunc[[1]](Xt = data, l = blocklen, r = bootrep, irreg = misstr) 
    if (interpolation == "Stineman"){
      z <- na.interpolation(boot.result, option = "stine")
      RegularSeries <- apply(z , 1, mean , na.rm = TRUE)
    } else if (interpolation == "Kalman"){
      z <- na.kalman(x = boot.result, model = "auto.arima", smooth = TRUE)
      RegularSeries <-apply(z , 1, mean , na.rm = TRUE)
    } else if (interpolation == "WMA"){
      z <- na.ma(x = boot.result, k = 4, weighting = "exponential")
      RegularSeries <-apply(z , 1, mean , na.rm = TRUE)
    }
  }else{    #combine with ivs
    
    boot.result = BootFunc[[1]](Xt = data, l = blocklen, r = bootrep, irreg = misstr)
    if (interpolation == "Stineman"){
      z <- na.interpolation(boot.result, option = "stine")
      RegularSeries <-apply(z , 1, mean , na.rm = TRUE)
    } else if (interpolation == "Kalman"){
      z <- na.kalman(x=boot.result, model = "auto.arima", smooth = TRUE)
      RegularSeries <- apply(z , 1, mean , na.rm = TRUE)
    } else if (interpolation == "WMA"){ 
      z <- na.ma(x=boot.result, k = 4, weighting = "exponential")
      RegularSeries <- apply(z , 1, mean , na.rm = TRUE)
    }
  }
  imputedSeries[listNA] <- RegularSeries[listNA]
  return(
    imputedSeries
  )
}