library(bio3d)
library(dplyr)

soglia <- 1000
flag_catena <- 'h'
setwd("Desktop/Università/Tesi/")
path.proteine <-"Python/Dati/proteine.csv"
proteine <- read.csv(path.proteine)
n.row <- nrow(proteine) #numero di righe
url <- "https://ftp.rcsb.org/pub/pdb/data/structures/divided/pdb"

CHAIN <-  "A"
contatore <- 0
for (i in 1:n.row){
  print(i)
  
  classe <- proteine[i,3]
  classe <- strsplit(classe, ".", fixed = TRUE)
  classe <- classe[[1]][1]
   
  idcentrale <- proteine[i,2]
  id <- proteine[i,1]
  code <- paste("pdb", id, ".ent.gz", sep = "")
  #Creazione URL
  url.finale <- paste(url, idcentrale, code, sep = "/")
  
  nome <- paste(id,"csv", sep = ".")
  #Creazione percorso
  path <- paste("R/Rdati/Matrici Adiacenza_otto_classi",nome, sep = '/')
  #print(path)
  #SOLO TEMPORANEO (IF con controllo sull'id)
  
  if(!file.exists(path) 
     && id != "1c0h" && id != "3n48"  && id !="2dfq" && id != "4i0v" 
     && id != "2dfo" && id != "2hbi" && id != "1q9q" && id != "1q9r" 
     && id != "1q9v" && id != "1q9t" && id != "2x0m" && id != "2yc6"
     && id != "2yc6" && id != "2yc7" && id != "2yc8" && id != "2ck9"
     && id != "1l7w" && id != "1nt7" && id != "6gwn" && id != "2v6d"
     && id != "3kgo" && id != "3kgn" && id != "3kgm" && id != "1l7u"
     && id != "2qkg" && id != "4lxe" && id != "3wmo" && id != "3wmn"
     && id != "2uwg" && id != "2uwz" && id != "1bxj" && id != "1lqr"
     && classe == flag_catena && contatore <= soglia){
    contatore <- contatore + 1
    print(contatore)
    
    print("Scarico le informazioni")
    
    pdb <- read.pdb(url.finale,ATOM.only = T) 
    
    #Select only ATOM records
    inds <- atom.select(pdb,type="ATOM",chain=CHAIN)
    
    #Get contact map matrix
    contact.map <- cmap(pdb, dcut=10, scut=0, inds=inds, mask.lower=F)
    contact.map[is.na(contact.map)] <- 0
    contact.map <- apply(contact.map,2,as.character)
    
    #Get the list of aminoacids
    list.amino <- pdb$atom[inds$atom,c("resno","resid")]
    by_resno <- list.amino %>% group_by(resno)
    res.table <- by_resno %>% summarise(res=unique(resid))
    list.amino <- res.table$res
    
    #Annotate rows and columns of contact map with aminoacid names
    contact.map <- rbind(list.amino,contact.map)
    #contact.map <- cbind(c("Residue",list.amino),contact.map)
    #View(contact.map)
    
    
    print("Salvataggio Matrice adiacenza")
    write.csv2(contact.map, file = path, quote=FALSE, row.names=T)
    
  }
}

