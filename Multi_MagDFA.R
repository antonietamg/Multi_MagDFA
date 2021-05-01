##################################
##################################
##################################
# MULTI MagDFA FUNCTION
##################################
##################################
##################################

path = "my_path"  # path to the folder with Multi_DFA.R file
setwd(path)

mdfa_function <- function(x, a_val, pol, sec){
  
  mdfa_function_final <- function(x){
    x_series <- data.frame(x[ , 2 : ncol(x)])
    t_series <- as.numeric(x[ , 1])
    integrate.serie <- function(x){     
      p <- cumsum((x) - mean(x))
      return(p)
  }
    
  integrate.serie <- cmpfun(integrate.serie)
  serie_int <- lapply(x_series, function(x) integrate.serie(x))
    
  variables.dfa <- function(x, y){
    tam = nrow(x)  
    
    if(is.null(tam) == TRUE){
        print("nrow to length")
        tam =length(x)
    }
    
    n = 4  
    tmax = trunc((tam / n), digits = 0)
    tmin = 10
    si = (tmax - tmin + 1)
    fu = (tmax - tmin + 1)
    N = tam
    nu=0
    dats <- sec
    ndats <- length(sec)
    dats <- as.numeric(sec)
      
    for (s in 1:ndats){
        nu[s] = (trunc(tam / dats[s]) - 1) #se calcula cuantos bloques de s elementos pueden obtenerse.
        Fn = tam
    }
    
    my_list <- list("tmax"=tmax, "tmin"=tmin, "N"=N, "dats"=dats, "ndats"=ndats, "Fn"=Fn, "si"=si, "fu"=fu, "nu"=nu, "s"=s )
    return(my_list)
  }
  
  variables.dfa <- cmpfun(variables.dfa)
  variables <- variables.dfa(x = x_series, y = t_series)
    
  blockf  <- function(w, x, s = s, k = k){
    block = x[(k * w[s] + 1) : (k * w[s] + w[s])]
  }
  
  blockf <- cmpfun(blockf)
  
  blockft  <- function(w, y, s = s, k = k){
    blockt = y[(k * w[s] + 1) : (k * w[s] + w[s])]
  }
  
  blockft <- cmpfun(blockft)
    
  dfa_fluctuation <- function(x){  
    Fs <- list()
    Fss =numeric()
      for (j in 1 : length(variables[["dats"]])){
        q <- c(1 : length(variables[["dats"]]))
        nu = variables[["nu"]]
        knu = nu[j]
        y = t_series
        z = knu
        w = variables[["dats"]]
        s = q[j]
        dfb <- list()
        dfbt <- list()
          for (k in 0:z){
            bl <- blockf(w, x, s = s, k = k)
            nam = sprintf("B%s_%s_%s", 1, w[s], k)
            dfb[[nam]] <- bl
          }
          for (k in 0:z){
            blt <- blockf(w, y, s = s, k = k)
            nam = sprintf("B%s_%s_%s", 1, w[s], k)
            dfbt[[nam]] <- blt
          }
        
        coef1 <- mapply(x = dfbt, y = dfb, FUN = function(x, y) polyfit(x, y, pol)) # el numero es el grado del polinomio
        coef1 <- data.frame(coef1)
        coef1 <- as.list(coef1)
        evalp <- mapply(x = coef1, y = dfbt, FUN = function(x, y) polyval(x, y))
        evalp <- data.frame(evalp)
        evalp <- as.list(evalp)
        k_var <- mapply(x = dfb, y = evalp, FUN = function(x, y) mean((x-y)^2))
        yevalp <- as.vector(k_var)
        Fnx <- yevalp
        Fss[[j]] <- sum(Fnx)
        Fs <- Fss
      }
      return(Fs)
    }
    
    dfa_fluctuation <- cmpfun(dfa_fluctuation)
    funcfn <- lapply(serie_int, function(x) dfa_fluctuation(x))
    
    sifu.alpha.dfa <- function(x, v, vari, l){
      lon = length(x)
      funcntot = 0
        for(i in 1:lon){
          serie_f <- x[[i]]
          funcntot <- funcntot+ sapply(serie_f, as.numeric)
        }
      v=v-1
      funcntot <- funcntot / (vari[["nu"]])
      Fnn = (sqrt(funcntot)) / (vari[["N"]])
      t1 <- vari[["dats"]][1]
      tt <- tail(vari[["dats"]], n = 1)
      si = vari[["si"]]
      fu = vari[["fu"]]
        for (s in 1 : length(sec)){
          si[s] = log10(sec[s])
          fu = log10(Fnn)
        }
      A = polyfit(si, fu)  #alpha value
      my_list <- list("A" = A, "si" = si, "fu" = fu)
      return(my_list)
    }
    l = nrow(x)
    sifu.alpha.dfa <- cmpfun(sifu.alpha.dfa)
    sifu <- sifu.alpha.dfa(x = funcfn, v = ncol(x), vari = variables, l)
    
    bin_func <- function(x){
      alpha_val <<- x[[1]]
      sx = x[[2]]
      fx = x[[3]]
      dplot <- data.frame()
      dplot <- data.frame(cbind(sx, fx))
      return(dplot)
    }
    bin_func <- cmpfun(bin_func)
    dplot_serie <- bin_func(sifu)
    return(dplot_serie)
  }
  
  mdfa_function_final <- cmpfun(mdfa_function_final)
  
  
  ###### INCREMENT SERIES
  
  l = nrow(x)
  nn = 0
  
  # IF ALPHA_VAL < 0.5 
  # "For anticorrelated series, the original series Xi may be treated as the increment series..."
  # Ashkenazy(2002), p.3
  
  if(a_val > 0.5){
    print("generating increment series...")
    serieinc <- function(x){
      serie_inc <- data.frame(2 : nrow(x))
      seriei = numeric()
      ttt = ncol(x)
        for(i in 2 : ttt){
          serx = sapply(x[,i], as.numeric)
          seriei <- diff(serx, 1)
          serie_inc <- cbind(serie_inc, seriei)
        }
      return(serie_inc)
    }
    serieinc <- cmpfun(serieinc)
    serieinc1 <- serieinc(x)
    x <- serieinc1
    nn = 1
  }
  
  # "In the case that the increment series is correlated with DFA exponent a > 0.5 it is neccessary to differenciate
  # the series until it becomes anticorrelated with exponent a < 0.5". Ashkenazy(2002), p.3 
  
  if(nn == 1){
    source("Multi_DFA.R", echo = TRUE)
    dfa_function(x, pol = pol, sec = sec) 
    av = alpha_val[1]
    print(alpha_val[1])
    
  while (av > 0.5){  #check alpha_value
    serieinc1 <- serieinc(x)
    x <- serieinc1
    dfa_function(x, pol = pol, sec = sec)
    av = alpha_val[1]
    print(alpha_val[1])
    print("alpha > 0.5")
    }
  }
  ###### 
  
  ###### MAGNITUDE SERIES
  
  integrate.serie <- function(x) {       #funcion para integrar la serie
    p <- cumsum((x) - mean(x))
    return(p)
  }
  
  integrate.serie <- cmpfun(integrate.serie)
  
  seriesms <- function(x){ 
    mag_s <- abs(x)
    mag_sf <- data.frame(1 : nrow(mag_s))
      for(i in 2 : ncol(mag_s)){
        mag_sx <- sapply(mag_s[ ,i], as.numeric)
        mag_sxx <- mag_sx - mean(mag_sx)
        mag_sf <- cbind(mag_sf, mag_sxx)
      }
    seriesms <- cmpfun(seriesms)
    mag_s2 <- data.frame(1 : nrow(mag_s))
    series = numeric()
    seriem = numeric()
    ttt = ncol(x)
      for(i in 2 : ttt){
        seriey = sapply(mag_sf[ ,i], as.numeric)
        serieyy <- integrate.serie(seriey)
        mag_s2 <- cbind(mag_s2, serieyy)
      }
    
    my_list <- list("mag_s2"=mag_s2)
    return(my_list)
  }
  series_mag.sig <- seriesms(x)
  ddd <- series_mag.sig[["mag_s2"]] 
  ddd <- mdfa_function_final(ddd)
  return(ddd)
}
mdfa_function <- cmpfun(mdfa_function)

##################################
##################################
##################################
##################################
# END MagDFA FUNCTION
##################################
##################################
##################################