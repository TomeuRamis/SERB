
#Dado una red bayesiana con conocimiento experto quiero tomar la decisión
#Esta decisión se obtine a partir de la "predición" que devuelve la red bayesiana dada cierta evidencia

#Objectivo: cargar red bayesiana, crear evidencias y obtener la respuesta de un nodo

#install.packages("bnlearn")
#library(bnlearn)
#library(tidyverse) #necesaria para ordenación
#
# #Cargamos la BN Asia con nombre bn
# load("R/asia.rda")
# load("R/sachs.rda")
# load("R/water.rda")
# #Mostramos la red cargada
# #install.packages("graphviz")
# #library(graphviz)
# graphviz.plot(bn)

#Function to generate the KB of a BN
extraer_niveles <- function(modelo, nod){
  err = TRUE
  niveles <- list() #Lista con los niveles de cada nodo
  for(i in 1:length(modelo)){ #Para cada uno de los nodos de la BN
    if(names(modelo[i]) != nod){ #Si se trata del nodo a predeci lo saltamos
      dimn <- dimnames(modelo[[i]]$prob)
      if(length(dimn) > 1){
        niveles[[length(niveles)+1]] <- dimn[[1]]
      }else{
        niveles[[length(niveles)+1]] <- dimn[[1]]
      }
    }else{
      err = FALSE
    }
  }

  #Error en caso de que el nodo indicado no se encuentre
  if(err){
    stop("Error during KB creation: the node specified was not found in the model.")
  }
  #print(paste("Niveles obtenidos: ",niveles))
  return(niveles)
}

generar_evidencia <- function(bn, node){
  n <- extraer_niveles(bn, node) #Obtenemos los niveles de cada nodo

  kb <- do.call(expand.grid, n) #Generamos todas las posibles combinaciones
  #Añadimos los nombres a las columnas
  names <- bnlearn::nodes(bn)
  colnames(kb) <- names[!names == node]
  kb_o <- dplyr::arrange_all(kb)  #Reordenamos la tabla
  #print(kb_o)
  return(kb_o)
}

#' @export
createKB <- function(BN, node){

  ev <- generar_evidencia(BN, node) #Generamos la evidencia

  #Creamos vector de salida
  D <- data.frame(matrix(ncol = 2, nrow = nrow(ev)))
  colnames(D) <- c("Decision", "Prob")

  #Iteramos por cada fila de tabla de conocimiento
  for( i in 1:nrow(ev)){
    P <- predict(BN, node = node, data = ev[i,], method = "bayes-lw", prob=TRUE)
    D[i,"Decision"] <- P
    D[i, "Prob"] <- max(attr(P,"prob"))
  }

  #Combinamos la evidencia con la decisión para obtener la Base de conocimiento
  KB <- cbind(ev,D )

  #Deifinimos Decision como factor y le asignamos los valores del modelo
  KB$Decision <- as.factor(KB$Decision)
  levels(KB$Decision) <- dimnames(BN[[node]]$prob)[[1]]

  #Eliminamos NA
  #KB <- na.omit(KB)
  #Eliminamos la columna predecida
  #KB <- subset(KB, select=-c(node))
  return(KB)
}

#' @export
writeKB <- function(file, KB){
  #Abrimos fichero de output y volcamos la tabla
  file.create(file)
  write.table(KB, file, sep=",", col.names = TRUE, row.names = FALSE)
  }

#' @export
readKB <- function(file){
  KB <- read.table(file, sep=",", header=TRUE)
  #Separate and drop last column, should allways be prob
  prob <- KB[ , ncol(KB)]
  rKB <- KB[1:(ncol(KB)-1)]
  #Aplicamos factores a todas las columnas (menos a prob)
  rKB[colnames(rKB)] <- lapply(rKB[colnames(rKB)], factor)
  #Juntamos de nuevo prob
  KB <- cbind(rKB,prob )
  return(KB)
}
#
# KB <- createKB(bn, "lung")
# #KB <- createKB(bn, "CBODD_12_45")
# writeKB("output_asia.txt",KB)
#
# n_KB <- readKB("output_asia.txt")
#
# #predict(bn, node = "tub", data = ev[1,], method = "bayes-lw", prob=TRUE)
# #predict(bn, node = "tub", data = evidence[1,], method = "bayes-lw", prob=TRUE)

