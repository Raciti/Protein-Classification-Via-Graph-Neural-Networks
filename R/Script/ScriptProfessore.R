library(bio3d)
pdb <- read.pdb( system.file("examples/hivp.pdb", package="bio3d") )
#Get the list of aminoacids
list.amino <- pdb$atom$resid
#Get contact map matrix
ref.cont <- cmap(pdb, dcut=10, scut=0) # cercare la catena 

#Modifiche fatte da me 
#Salvataggio della matrice

write.csv2(ref.cont, file = "Desktop/Università/Tesi/Rdati/datiR.csv", quote=FALSE, row.names=T)

ref.cont


#######
pdb <- read.pdb("Downloads/pdb2w9l.ent")
#Get the list of aminoacids
list.amino <- pdb$atom$resid
list.chain <- pdb$atom$chain
for(i in 1:length(list.chain)){
  if(list.chain[i] != "A"){
    index <- i
    break
  }
}
print(index)
#Get contact map matrix
proteina <- read.pdb("../../../../../Downloads/pdb100d.ent")
list.amino <- proteina$atom$resid
inds <- atom.select(proteina, "protein", chain="A")
ref.cont <- cmap(proteina$xyz[inds$xyz], dcut=10, scut=0)
write.csv2(ref.cont, file = "Desktop/Università/Tesi/Rdati/Matrici Adiacenza/2hb1.csv", quote=FALSE, row.names=T)


