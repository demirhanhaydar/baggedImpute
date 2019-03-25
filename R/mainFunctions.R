#' @import stats
#' 

GenMBB <- function(Xt, l = 2, r = 100, irreg = c(NULL, "ams", "ivs")){
  
  if(is.null(Xt)){
    stop("Error: data need to be included!")
  }
  
  if(is.null(irreg) | irreg == "" | irreg == "ams"){
    n <- length(Xt)
    if(l>n){
      stop("Error: block length must be less than or equal to the series length.")
      
    }
    
    numBlocks <- ceiling(n-l+1)
    Xt_matrix <- matrix(nrow = l, ncol = numBlocks)
    
    for(i in 1:numBlocks){
      for(j in 1:l){
        Xt_matrix[j,i] = Xt[(i-1)+j]
      }
    }

    ss <- ceiling(n/l)
    hxt_mat <- c(1:numBlocks) #--> to help us make the index to be resampled
    
    return(
      ts( replicate(r,as.vector(Xt_matrix[,sample(hxt_mat, size = ss, replace = TRUE)])) )
    )
    
  }else if(irreg == "ivs"){
    
    n <- length(na.exclude(Xt))
    
    if(l>n){
      stop("Error: block length must be less than or equal to the series length.")
    }
    
    numBlocks = length(which(!is.na(Xt))) - l + 1 
    varlist = list()
    
    for(i in 1:numBlocks){
      j <- 0
      k <- 1
   
      Xt_first_obs_position <- which(!is.na(Xt))[i]
      
      varlist[[i]] <- Xt[Xt_first_obs_position]
      
      while(j < (l - 1)){
        if(is.na(Xt[Xt_first_obs_position+k])){
          varlist[[i]][k+1] <- NA
        }else{
          varlist[[i]][k+1] <- Xt[Xt_first_obs_position+k]
          j <- j+1
        }
        k <- k+1
      }
    }
    
    bootList <- list()
    
    for(i in 1:r){
      bootList[[i]] <- unlist(varlist[sample(x = c(1:numBlocks), size = 1, replace = TRUE)])
      while(length(bootList[[i]])<length(Xt)){
        bootList[[i]] <- append(x = bootList[[i]], values = unlist(varlist[sample(x = c(1:numBlocks), size = 1, replace = TRUE)]), after = length(bootList[i]))
      }
    }
    return(
      ts(matrix(unlist(bootList), nrow = length(Xt), ncol = r))
    )
    
  }else{
    stop("Please, correctly specify irreg input with one of NULL, ams, or ivs.")
  }
}


GenNBB <- function(Xt, l = 2, r = 100, irreg = c(NULL, "ams", "ivs")){
  #Read the data, and store the data length
  if(is.null(Xt)){
    stop("Error: data need to be included!")
  }
  
  if(is.null(irreg) | irreg == "" | irreg == "ams"){

    n <- length(Xt)
    if(l>n){
      stop("Error: block length must be less than or equal to the series length!")
      
    }

    numBlocks <- floor(n/l)
    
    #Generating matrix data to store the block values 
    Xt_matrix <- matrix(nrow = l, ncol = numBlocks)
    
    for(i in 1:numBlocks){
      for(j in 1:l){
        Xt_matrix[j,i] <- Xt[(l*(i-1))+j]
      }
    }
    
    ss <- ceiling(n/l)
    hxt_mat <- c(1:numBlocks) 
    return(
      ts(replicate(r,as.vector(Xt_matrix[,sample(hxt_mat, size = ss, replace = TRUE)])))
    )
    
  }else if(irreg == "ivs"){
    
    n <- length(na.exclude(Xt))
    if(l>n){
      stop("Error: block length must be less than or equal to the series length!")
    }
    
    numBlocks = floor(n/l) 

    varlist = list()
    for(i in 1:numBlocks){
      j = 0
      k = 1
      
      Xt_first_obs_position = which(!is.na(Xt))[l*(i-1)+1]
      
      varlist[[i]] <- Xt[Xt_first_obs_position]
      
      while(j < (l - 1)){
        if(is.na(Xt[Xt_first_obs_position+k])){
          varlist[[i]][k+1] <- NA
        }else{
          varlist[[i]][k+1] <- Xt[Xt_first_obs_position+k]
          j <- j+1
        }
        k <- k+1
      }
    }
    bootList <- list()
    
    for(i in 1:r){
      bootList[[i]] <- unlist(varlist[sample(x = c(1:numBlocks), size = 1, replace = TRUE)])
      while(length(bootList[[i]])<length(Xt)){
        bootList[[i]] <- append(x = bootList[[i]], values = unlist(varlist[sample(x = c(1:numBlocks), size = 1, replace = TRUE)]), after = length(bootList[i]))
      }
    }
    return(
      ts(matrix(unlist(bootList), nrow = length(Xt), ncol = r))
    )
    
  }else{
    stop("Please, correctly specify irreg input with one of NULL, ams, or ivs.")
  }
  
}

GenCBB <- function(Xt, l = 2, r = 100, irreg = c(NULL, "ams", "ivs")){
  
  #Read the data, and store the data length
  if(is.null(Xt)){
    stop("Error: data need to be included")
  }
  
  if(is.null(irreg) | irreg == "" | irreg == "ams"){
    n <- length(Xt)
    if(l>n){
      stop("Error: block length must be less than or equal to the series length!")
    }
    Xt_matrix <- matrix(nrow = l, ncol = n)
    Xtr <- c(Xt, Xt[1:(l-1)])
    
    for(i in 1:n){
      for(j in 1:l){
        Xt_matrix[j,i] = Xtr[(i-1)+j]
      }
    }

    ss <- ceiling(n/l)
    
    hxt_mat <- c(1:n) 
    return(
      ts(replicate(r,as.vector(Xt_matrix[,sample(hxt_mat, size = ss, replace = TRUE)])))
    )
    
  }else if(irreg == "ivs"){
    
    n <- length(na.exclude(Xt))
    
    if(l>n){
      stop("Error: block length must be less than or equal to the series length!")
    }
    
    varlist = list()
    
    Xtr <- c(Xt, Xt[1:(which(!is.na(Xt))[l-1])])
    
    for(i in 1:n){
      j <- 0
      k <- 1
      
      Xtr_first_obs_position = which(!is.na(Xtr))[i]
      
      varlist[[i]] <- Xtr[Xtr_first_obs_position]
      
      while(j < (l - 1)){
        if(is.na(Xtr[Xtr_first_obs_position+k])){
          varlist[[i]][k+1] <- NA
        }else{
          varlist[[i]][k+1] <- Xtr[Xtr_first_obs_position+k]
          j <- j+1
        }
        k <- k+1
      }
    }
    
    bootList <- list()
    
    for(i in 1:r){
      bootList[[i]] <- unlist(varlist[sample(x = c(1:n), size = 1, replace = TRUE)])
      while(length(bootList[[i]])<length(Xt)){
        bootList[[i]] <- append(x = bootList[[i]], values = unlist(varlist[sample(x = c(1:n), size = 1, replace = TRUE)]), after = length(bootList[i]))
      }
    }
    return(
      ts(matrix(unlist(bootList), nrow = length(Xt), ncol = r))
    )
    
  }else{
    stop("Please, correctly specify irreg input with one of NULL, ams, or ivs.")
  }
}

