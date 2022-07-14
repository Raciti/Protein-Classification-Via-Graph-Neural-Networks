library(bio3d)

path.proteine <-  "Desktop/Università/Tesi/Python/Dati/proteine.csv"
proteine <- read.csv(path.proteine)
n.row <- nrow(proteine) #numero di righe
url <- "https://ftp.rcsb.org/pub/pdb/data/structures/divided/pdb"

for (i in 1:n.row){
  print(i)

  idcentrale <- proteine[i,2]
  id <- proteine[i,1]
  code <- paste("pdb", id, ".ent.gz", sep = "")
  #Creazione URL
  url.finale <- paste(url, idcentrale, code, sep = "/")
  
  nome <- paste(id,"csv", sep = ".")
  #Creazione percorso
  path <- paste("Desktop/Università/Tesi/R/Rdati/Matrici Adiacenza",nome, sep = '/')
  #print(path)
  #SOLO TEMPORANEO (IF con controllo sull'id)
  if(id != "4i0v"  && !file.exists(path) && id != "1c0h" && id != "3n48"  && id !="2dfq"
     && id !="2dfo" && id != "2hbi"){
    print("Scarico le informazioni")
    proteina <- read.pdb(url.finale)
    
    #print("Lista amino acidi")
    #print("Creazione Matrice adiacenza")
    
    inds <- atom.select(proteina, "protein", chain="A") #Selezione degli aminoacidi della catena A
    
    
    ref.cont <- cmap(proteina$xyz[inds$xyz], dcut=10, scut=0, mask.lower = FALSE) #Matrice Adiacenza, mask.lower non mette gli NA
    #ref.cont <- cmap(proteina, dcut=10, scut=0)
    
    #Eliminazione self-loop
    diag(ref.cont) <- 0
  
    list.amino <- proteina$atom$resid         #Lista aminoacidi
    list.amino <-list.amino[1:nrow(ref.cont)] #Lista aminoacidi ridimensionata in base a quelli selezionati
    
    
    ref.cont <- rbind(list.amino, ref.cont)
    
    
    print("Salvataggio Matrice adiacenza")
    write.csv2(ref.cont, file = path, quote=FALSE, row.names=T)
    break
  }
}

