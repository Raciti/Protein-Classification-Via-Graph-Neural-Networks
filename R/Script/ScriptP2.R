library(bio3d)
library(dplyr)

setwd("Downloads/")
#pdb1iv8.ent
PDB.ID <- "pdb1iv8"
CHAIN <- "A"

#Read PDB 3D structure
pdb <- read.pdb(paste0(PDB.ID,".ent"),ATOM.only = T)

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
contact.map <- cbind(c("Residue",list.amino),contact.map)
#View(contact.map)

#Write contact map
write.table(contact.map,paste0(PDB.ID,".txt"),quote=F,sep="\t",row.names=F,col.names=F)
