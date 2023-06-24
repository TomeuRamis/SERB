#Dada una base de conocimiento queremos aplicar un algoritmo de optimización que reduzca el número de count.rules al máximo (por ahora)
#library(tidyverse) #necesaria para ordenación

#Hamming distance function
H <- function(x, y){
  sum(x != y)
}

#Number of attributes between the initial (X) and the final (Y) position of a permuted attribute i
Ri <- function(x, y, i){
  #obtenemos nombre de la columna i en x
  name = x[i]
  #obtenemos posición de la columna i en y
  j = grep(name, y)
  #Calculamos atributos intermedios
  return(i - j - 1)
}

#G Distance function between 2 basesX and Y
G <- function(x, y){
  #G(B,B)=0
  if(identical(x,y)){
    return(0)
  }

  v <- as.vector(c())
  sum = 0

  for(i in length(x):1){ #For each attribute from right to left
    #print(paste("Columna ",x[i], " y columna ", y[i]))
    if(x[i] != y[i]){ #If it has been permuted
      #print("No coinciden!")
      firstDiff = i-1
      if(x[i] %in% v){ #If the attribute has been counted ignore
        #print(paste("La columna ",x[i], " ya había sido contada"))
        next
      }
      #print(paste("Vector: ",v))
      #Add the attributes to the treated vector
      v <- append(v, x[i])
      v <- append(v, y[i])
      #Add the number of atributtes between the permutated columns
      sum = sum + Ri(x, y, i)
      #print(paste("Ri of ",x[i]," and ",y[i]," at ", i, " is ", Ri(x, y, i)))
    }
  }
  #print(paste("fistDiff: ",firstDiff, ", sum: ", sum,", Hamming:", H(x,y)))
  return(firstDiff + sum + H(x,y) - 1)
}

#Permutación de dos columnas
swtch <- function(x,i,j) {x[c(i,j)] <- x[c(j,i)]; x}

#Intercambiar dos columnas y reordenar los datos
SwapColumns <- function(kb, x, y){
  if(x == y){
    return(0) #ERROR
  }
  oC <- colnames(kb)

  nC <- swtch(oC, x, y)
  newKb <- kb[nC]
  #Recoger info de ordenación
  aux <- newKb[,1:(length(newKb)-2)]
  #print(aux)
  #Concatenar valores en forma de string
  f <- apply(aux, 1, paste, collapse = "_")
  #Ordenar por orden alfabetico basandote en f
  orderedKb <- newKb[order(f, decreasing = FALSE),]
  #print(orderedKB)
  return(orderedKb)
}

# This function returns TRUE wherever elements are the same, including NA's,
# and FALSE everywhere else.
compareNA <- function(v1,v2) {
  same <- (v1 == v2) | (is.na(v1) & is.na(v2))
  same[is.na(same)] <- FALSE
  return(same)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#A drop-in replacement for {base::rle()} that treats all NAs as identical
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rle2 <- function (x)  {
  stopifnot("'x' must be a vector of an atomic type" = is.atomic(x))

  n <- length(x)
  if (n == 0L) {
    return(structure(list(
      lengths = integer(), values = x)
    ), class = 'rle')
  }

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Where does next value not equal current value?
  # i.e. y will be TRUE at the last index before a change
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  y <- (x[-1L] != x[-n])

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Since NAs are never equal to anything, NAs in 'x' will lead to NAs in 'y'.
  # These current NAs in 'y' tell use nothing - Set all NAs in y to FALSE
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  y[is.na(y)] <- FALSE

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # When a value in x is NA and its successor is not (or vice versa) then that
  # should also count as a value change and location in 'y' should be set to TRUE
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  y <- y | xor(is.na(x[-1L]), is.na(x[-n]))

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Any TRUE locations in 'y' are the end of a run, where the next value
  # changes to something different
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  i <- c(which(y), n)

  structure(list(
    lengths = diff(c(0L, i)),
    values  = x[i]
  ), class = 'rle')
}

#Contar la cantidad de count.rules de una base de conocimiento. Si v = TRUE (verbose), sacamos por consola las count.rules.
#' @export
count_rules <- function(KB, verbose){
  #rle$values devuelve una lista con cada conjunto de decisiones iguales contiguas
  rle <- rle2(KB$Decision)
  if(verbose){
    s = ""
    sum = 0
    for(i in rle$lengths){
      s = paste(s, "< (",KB[i+sum,1] )
      for(j in 2:(length(KB[i+sum,])-2)){
        s = paste(s,",", KB[(i+sum),j])
      }
      s = paste(s, ")", KB$Decision[(i+sum)], "|")
      sum = sum + i
    }
     print(s)
  }
  return(length(rle$values))
}

#Búsqueda tabú tomando como herística 'count.rules'
#' @export
tabu <- function(KB, m, maxtabu){
  #Inicializamos
  end = FALSE #condición de terminación
  mp = 0
  best <- KB
  bCandidato <- KB
  p = length(KB)-2 #permutaciones posibles
  t <- vector() #tabú
  r <- matrix(nrow = p, ncol=p)

  while(!end){
    #Generamos vecinos
    for(i in 1:(p-1)){
      for(j in (i+1):p){
        r[i,j] <- G(colnames(bCandidato),swtch(colnames(bCandidato),i,j))
        #r[i,j] <- count.rules(SwapColumns(best, i, j))
      }
    }

    #Buscar mejor vecino que no sea tabú
    n <- which(r == min(r, na.rm = TRUE), arr.ind = TRUE)
    base <- swtch(best, n[1,1],n[1,2])
    b <- paste(base, collapse = "_")
    while(b %in% t ){
      r[n[1,1],n[1,2]] <- NA
      n <- which(r == min(r, na.rm = TRUE), arr.ind = TRUE)
      base <- swtch(best, n[1,1],n[1,2])
      b <- paste(base, collapse = "_")
    }

    if(!(b %in% t)){ #si no está en el vector tabú
      bCandidato <- SwapColumns(bCandidato, n[1,1], n[1,2])
      if(count_rules(bCandidato, FALSE) < count_rules(best, FALSE)){ #Si es mejor que el actual mejor
        best <- bCandidato
        #print(bCandidato)
      }

      t <- append(t, b)
      if(length(t) > maxtabu){ #If tabú list is at max capacity remove first element
        t <- t[-1]
      }
    }

    #Condición de fin
    if(mp >= m){
      end = TRUE
    }else{
      mp = mp + 1
    }
  }
  return(best)
}

#Simulated Annealing
#' @export
sima <- function(KB, T0, k, c){
  #Inicializamos
  end = FALSE
  Temp = T0
  Tf = 0.1
  p = length(KB)-2 #permutaciones posibles
  r <- matrix(nrow = p, ncol=p)

  best <- KB

  while(!end){
    #Generamos vecinos
    for(i in 1:(p-1)){
      for(j in (i+1):p){
        r[i,j] <- G(colnames(best),swtch(colnames(best),i,j))
      }
    }

    #Buscar mejor vecino
    n <- f(T0, Tf, Temp, r)
    #print(paste("Distancia G sugerida: ", n))
    N <- which(r == n, arr.ind = TRUE)
    while(all(is.na(N))){
      np <- ceiling(n)
      if(np == n){
        np <- (n + 1)
      }
      n <- np
      N <- which(r == n, arr.ind = TRUE)
    }
    #print(paste("Distancia G elegida: ", n))
    bCandidato <- SwapColumns(best, N[1,1], N[1,2])

    #Es mejor que la solución actual?
    if(count_rules(bCandidato, FALSE) < count_rules(best, FALSE)){
      best <- bCandidato
    }else{
      ran <- runif(1) #Random number
      dif <- count_rules(bCandidato, FALSE) - count_rules(best, FALSE)-0.1
      #print(paste("Diferencia: ", dif))
      #print(paste("Chance = ", ran, ">",exp(-(dif/(k*Temp)))))

      if(ran > exp(-(dif/(k*Temp)))){ #Probabilidad aleatoria de hacer un "mal" movimiento
      #if(ran > exp(-(k*Temp))+0.25){
        best <- bCandidato
      }
    }
    Temp = Temp * c
    #print(paste("Temperatura actual: ", Temp))
    if(Temp < Tf){
      end = TRUE
    }
  }
  return(best)
}

f <- function(Tmax, Tmin, Tc, G){
  cG <- c(G)
  cG <- cG[!is.na(cG)]
  return(((Tc - Tmin)/(Tmax - Tmin))*max(cG)) # % de finalización X distancia máxima de G
}

#Función de selección del AG
seleccion <- function(r, q){
  return(q*(1-q)^(r-1))
}

#Dada una base origen (O) y una destino (D), devuelve las permutaciones necesarias para transformar O en D
permutaciones <- function(O, D){
  n <- length(O)
  per <- data.frame()
    for(i in 1:n){
      #print(paste(O[i]," ", D[i]))
      if(O[i] != D[i]){
        p <- match(O[i],D)
        per <- rbind(per, c(i,p))
        Dp <- D
        Dp[i] <- D[p]
        Dp[p] <- D[i]
        D <- Dp
    }
  }
  per <- per[nrow(per):1,]
  return(per)
}

#función de Ranking de AG
rank <- function(P, B){
  #Eliminamos temporalmente Prob
  B <- B[1:(length(B)-1)]
  #Obtenemos la base de conocimineto reducida
  v <- count_rules(B, FALSE) #nº count_rules B
  #print(v)
  KBr <- data.frame(matrix(ncol = ncol(B), nrow = 0)) #dataset con 1 entrada por regla
  colnames(KBr) <- colnames(B)
  for(i in 2:nrow(B)){
    #print(B[i,])
    #Gestionamos de forma especial los NA
    if((is.na(B[i,"Decision"]) && !is.na(B[i-1,"Decision"])) || (!is.na(B[i,"Decision"]) && is.na(B[i-1,"Decision"]))){
      #print(paste("Nueva fila: ", B[i-1,]))
      KBr[nrow(KBr)+1,] <- B[i-1,]
    }else if(!is.na(B[i,"Decision"]) && !is.na(B[i-1,"Decision"])){
      if(B[i,"Decision"] != B[i-1,"Decision"]){
        #print(paste("Nueva fila: ", B[i-1,]))
        KBr[nrow(KBr)+1,] <- B[i-1,]
      }
    }
  }
  #print(KB)
  #Calculamos el valor de cada individuo de la población
  Pv <- vector()
  Bcol <- colnames(B)[-length(colnames(B))]
  #print("Base actual:")
  #print(Bcol)
  for(i in 1:length(P)){
    #print(paste("length de P:",length(P)))
    per <- permutaciones(Bcol, P[[i]])
    #print(per)
    KBPi <- KBr
    if(length(per) != 0){
      for(j in 1:nrow(per)){
        #print(paste(per[j,1], " y ", per[j,2]))
        KBPi <- SwapColumns(KBPi, per[j,1],per[j,2])
      }
    }
    Pv[i] <- count_rules(KBPi, FALSE)
  }
  #print(Pv)
  #Ordenamos los individuos por su valor
  PRanked <- list()
  PRv <- vector()
  for(i in 1:length(P)){
    j <- which.min(Pv)
    PRanked[[i]] <- P[[j]]
    PRv[i] <- Pv[j]
    Pv[j] <- NA
  }

  return(list(PRanked, PRv))
}

#función de voting crossover
crossover <- function(Bf, Bm, Lf, Lm){
  a <- vector()
  for(i in 1:(length(Bf))){ #Contamos desde 1 hasta N, no de 0 a N-1 como en el paper!!!!!!!!!
    j <- i - 1
    k <- match(Bf[i],Bm)-1
    ai <- (j*Lm + k*Lf) / (Lf * Lm)
    a[i] <- ai
    #print(paste("a",Bf[i],": ",a[i]))
  }

  child <- vector()
  for(i in 1:length(a)){
    j <- which.min(a)
    child[i] <- Bf[j]
    a[j] <- NA
  }

  return(child)
}

#función de mutación
mutation <- function(B, m){
  for(i in 1:length(B)){
    ran <- runif(1) #Random number
    #print(paste("mutación si: ", ran," <= ", m))
    if(ran <= m){
      j <- sample.int(length(B)-1,1)
      if(j > i){
        j <- j + 1
      }
      B <- swtch(B, i, j)
    }
  }
  return(B)
}

generate_next_gen <- function(P, Pv, np, q, m){
  selection <- vector()
  for(i in 1:length(P)){
    selection[i] <- seleccion(i, q) #Calculamos las probabilidades de selección de cada individuo
  }

  H <- list()
  i <- 1
  while(length(H) < np){
    #print(paste("length H:",length(H)))
    Bf <- NA
    Bm <- NA
    Bfi <- 0
    Bmi <- 0

    na <- is.na(Bf)
    while(na){
      #Seleccionamos un individuo al azar
      #print(length(P))
      Bfi <- sample(1:length(P), 1)
      #print(paste("BFi:", Bfi))
      ran <- runif(1)
      #print(paste("runif; ", ran, " < ", selection[Bfi]))
      if(ran < selection[Bfi]){ #Si aleatoriamente se selecciona
        Bf <- P[[Bfi]] #Seleccionamos padre
        na <- FALSE
      }
    }
    na <- is.na(Bm)
    while(na){
      #Seleccionamos un individuo al azar
      #print(length(P))
      Bmi <- sample(1:length(P), 1)
      #print(paste("Bmi:", Bmi))
      ran <- runif(1)
      #print(paste("runif; ", ran, " < ", selection[Bmi]))
      if(Bmi != Bfi){ #Si no es el mismo individuo
        if(ran < selection[Bmi]){ #Si aleatoriamente se selecciona
          Bm <- P[[Bmi]] #Seleccionamos madre
          na <- FALSE
          }
      }
    }

    C <- crossover(Bf, Bm, Pv[Bfi], Pv[Bmi])
    c <- mutation(C, m)
    #print(C)
    H[[i]] <- C
    i <- i + 1

  }
  return(H)
}

#función para la generación de población inicial
initial_p <- function(Bo, n, m){
  P <- list()
  P[[1]] <- Bo
  for(i in 2:n){
    P[[i]] <- mutation(Bo, m)
  }
  return(P)
}

#' @export
genetic <- function(KBo, np, m, q, t){
  #Eliminamos la columna probabilidad
  #KBo <- KBo[1: ncol(KBo)-1 ]

  #Obtenemos base original a partir de la KB original (sin Decision)
  #Bo <- colnames(KBo)[-length(colnames(KBo))]
  Bo <- colnames(KBo[1:(length(KBo)-2)])
  iter = 0
  changes = TRUE

  #Generamos primera población
  Po <- initial_p(Bo, np, 0.5)

  #Ordenamos los individuos
  Aux <- rank(Po, KBo)
  P <- Aux[[1]] #Poblaciónn actual
  #print("P0:")
  #print(P)
  Pv <- Aux[[2]] #Fintness de cada individuo

  H <- list()

  while(iter < t && changes){

    #Generamos siguiente generación
    H <- generate_next_gen(P, Pv, np-2, q, m)
    #Introducimos los 2 mejores individuos a la siguiente generación
    H[[np-1]] <- P[[1]]
    H[[np]] <- P[[2]]

    #print(paste("Generación ", iter))
    #print(H)

    #Actualizamos la KB con la "mejor" de la generación anterior
    per <- permutaciones(Bo, P[[1]])
    if(length(per) != 0){
      for(j in 1:nrow(per)){
        KBo <- SwapColumns(KBo, per[j,1],per[j,2])
      }
    }
    #Actualizamos la base (sin Decision ni Prob)
    Bo <- colnames(KBo[1:(length(KBo)-2)])

    #Ordenamos los individuos
    Aux <- rank(H, KBo)

    #Actualizamos la población actual
    P <- Aux[[1]] #Poblaciónn actual
    #print("Población ranked:")
    #print(P)
    Pv <- Aux[[2]] #Fintness de cada individuo

    #Next iteration
    iter <- iter + 1
    #print(paste("iteration: ",iter))
    #Falta calcular si han sucedido cambios
  }

  return(KBo)
}

#A partir de una tabla de decisión generamos la KBM2L asociada
#' @export
createKBM2L <- function(KB){
  #print(v)
  rle <- rle2(KB$Decision)
  m <- c(mean(KB$Prob[1:rle$lengths[1]]))
  for(i in 2:length(rle$lengths)){
    rle$lengths[i] = rle$lengths[i-1] + rle$lengths[i]
    m[i] <- mean(KB$Prob[(rle$lengths[i-1]+1):rle$lengths[i]])
  }
  print(m)
  KBM2L <- KB[rle$lengths,]
  KBM2L$Prob <- m
  print(KBM2L)
  return(KBM2L)
}

reconstruir_KB <- function(kbm2l){
  lvls <- list()
  #Reconstruimos la KB
  for(i in 1:(length(kbm2l)-2)){
    l <- list()
    l <- levels(kbm2l[,i])
    lvls[[i]] <- l
  }
  nKB <- do.call(expand.grid, lvls)
  nKB <- dplyr::arrange_all(nKB)  #Reordenamos la tabla
  colnames(nKB) <- colnames(kbm2l)[1:(length(kbm2l)-2)]
  #Añadimos la Decisión y prob
  D <- data.frame(matrix(ncol = 2, nrow = nrow(nKB)))
  colnames(D) <- c("Decision", "Prob")
  D$Decision <- factor(D$Decision, levels=levels(kbm2l$Decision))
  inf <- 1
  found <- TRUE
  for(i in 1:nrow(kbm2l)){ #Para cada una de las count.rules
    #print(kbm2l[i,])
    for(r in inf:nrow(nKB)){ #Recorremos la tabla
      #print(nKB[r,])
      found <- TRUE
      for(j in 1:length(nKB)){ #Comparamos columna a columna si se trata del item de la KBM2L
        if(!compareNA(nKB[r,j], kbm2l[i,j])){ #Si no coinciden pasamos
          found <- FALSE
          break
        }
      }
      if(found){
        #print(nKB[r,])
        D[inf:r,"Decision"] <- kbm2l[i,"Decision"]
        D[inf:r,"Prob"] <- kbm2l[i,"Prob"]
        #print(D)
        inf <- r +1
        r <- inf
        break
      }
    }
  }
  KB <- cbind(nKB,D)
  return(KB)
}

#Dada una KBM2L y un item devuelve la parte fija
#' @export
fixed <- function(kbm2l, evidence){
  KB <- reconstruir_KB(kbm2l)
  kb_r <- KB[1:(length(KB)-2)]
  i <- 0
  equal = TRUE
  for(j in 1:nrow(kb_r)){
    for(k in 1:length(kb_r)){
      if( !compareNA(kb_r[j,k], evidence[1,k])){ #No me gusta este método tan lento
        equal = FALSE
        break
      }
    }
    if(equal){
      i <- j
      break
    }
    equal <- TRUE
  }
  inf <- 1
  #buscamos inf
  if(i != 1){
    for(k in (i-1):1){
      #Si la decisión anterior es distinta a la de nuestro item
      if(!compareNA(KB[k,]$Decision, KB[i,]$Decision)){
        inf <- k + 1
        break
      }
    }
  }
  sup <- nrow(KB)
  if(i != sup){#buscamos sup
    for(k in (i+1):nrow(KB)){
      #Si la decisión siguiente es distinta a la de nuestro item
      if(!compareNA(KB[k,]$Decision, KB[i,]$Decision)){
        sup <- k - 1
        break
      }
    }
  }
  fix <- length(KB)
  #Buscamos las columnas que coinciden
  for(j in 1:length(KB)){
    if(!compareNA(KB[inf,j], KB[sup,j])){
      fix <- j-1
      break
    }
  }
  return( colnames(KB)[0:fix]) #Devolvemos los 'fix' primeros nombres
}

#' @export
query.list <- function(KBM2L, evidence){
  if(!is.vector(evidence)){
    evidence <- as.vector(evidence)
  }
  KB <- reconstruir_KB(KBM2L)
  p <- vector()

  for(i in 1:nrow(KB)){
    found <- TRUE
    for(j in 1:length(evidence)){
      if(KB[i,j] != evidence[j]){
        found <- FALSE
        break
      }
    }
    if(found){
      p <- append(p, i)
    }
  }
  return(KB[p,] )
}

#
# #Load KB
# KB <- cargar_KB("output_asia.txt")
#
# #Preparación Test
# KB1 <- KB[,colnames(KB)]
# KB2 <- KB[,c("asia","tub","bronc","lung","smoke","either","xray", "dysp","Decision")] #Permutated 'smoke' with 'bronc'
#
# #Test Ri
# Ri(KB1, KB2, 5)
#
# #Test function
# G(KB1,KB2)
#
# #Test cambiar columnas
# nKB <- SwapColumns(KB, 1, 5)
# #Distancia G de las dos bases de conocimiento
# G(KB, nKB)
#
# #Test count.rules
# count.rules(KB, FALSE)
#
# #Test tabú
# rT <- tabu(KB, 100, 10)
# count.rules(rT, TRUE)
#
# #Test SA
# rSA <- sima(KB, 1000, 1, 0.99)
# count.rules(rSA, TRUE)
#
# #Test crossover
# Fa <- c("2", "1", "0", "3", "5", "4")
# Mo <- c("3", "2", "1", "4", "5", "0")
# Lf <- 30
# Lm <- 45
#
# crossover(Fa, Mo, Lf, Lm)
#
# #Test permutaciones
# O <- c("C", "B", "A", "D", "F", "E")
# D <- c("D", "C", "B", "E", "F", "A")
# permutaciones(O,D)
#
# Op <- swtch(O, 4, 6)
# Op <- swtch(Op, 3, 6)
# Op <- swtch(Op, 2, 3)
# Op <- swtch(Op, 1, 2)
# Op
#
# #Test mutación
# A <- c("A","B","C","D","E")
# mutation(A, 0.05)
#
# #Test población inicial
# initial_p(A, 100, 0.5)
#
# #Test rank
# #generamos tabla con todas las posibles combinaciones
# KBt = expand.grid(c("yes", "no"), c("yes", "no"),c("yes", "no"), c("yes", "no"),c("yes", "no"))
# colnames(KBt) <- A
# KBt$Decision <- sample.int(2,32, replace=TRUE)
# po <- initial_p(A,10,0.5)
# x <- rank(po,KBt)
# Ranked <- x[[1]]
# Values <- x[[2]]
#
# #Test de la optimización de Rank
# for(i in 1:length(Ranked)){
#   per <- permutaciones(A, Ranked[[i]])
#   KBPi <- KBt
#   if(length(per) != 0){
#     for(j in 1:nrow(per)){
#       KBPi <- SwapColumns(KBPi, per[j,1],per[j,2])
#     }
#   }
#   print(paste("count.rules de ", i,": ", count.rules(KBPi,FALSE)))
# }
#
# H <- generate_next_gen(Ranked, Values, 8, 0.5, 0.05)
#
# #Test GA
# Bo <- colnames(KB)
# KBGA <- genetic(KB, 20, 0.05, 0.5, 50)
# count.rules(KBGA,TRUE)
#
# KBM2L <- createKBM2L(KBGA)
# levels(KBM2L$Decision)
#
# KBreconstruida <- reconstruir_KB(KBM2L)
#
# #Test identify fixed part
# item <- KBGA[2,]
# item <- item[1:(length(item)-2)]
# item
# item2 <- KBGA[1,]
# item2 <- item2[1:(length(item2)-2)]
# item2
# item3 <- KBGA[3,]
# item3 <- item3[1:(length(item3)-2)]
# item3
# item = item2
# as.character(item)==as.character(item3)
#
# item <- c("no","no","no", "no", "no")
# fixed(KBGA, item)
#
# query.list(KBM2L, item)


