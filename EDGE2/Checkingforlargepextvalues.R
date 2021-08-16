load("D:/Documents/ResearchProject/Results/MASTER_ALLEDLOSS_FINAL.Rdata")

groups <- unique(master.all$group)

ACT <- subset(master.all, master.all$group =="Act")

ED <- list()


ED[c(1:6)] <- NA
names(ED) <- groups

for (i in 1:length(groups)) {
  ED[[i]] <- subset(master.all, master.all$group == groups[i])
}

save(ED, file = "AllEDLoss_incomplete.Rdata")



### checking for values >>1 
load("D:/Documents/ResearchProject/EDGE2/Amphibians/6000_pext.Rdata")
load("D:/Documents/ResearchProject/EDGE2/Mammals/7000_pext.Rdata")

for (i in 1:length(res.list)) {
  mac <- c()
  df <- res.list[[i]]
  
  for (j in 2:1001) {
   mac <- c(mac, max(df[1,j], na.rm = T))
  }
  
  if ( max(mac, na.rm = T) > 1) {
    mac <- subset(mac, mac >1)
    print(i)
    print(length(mac))
  }

}

# amphibians - 2,3,4 are too large
# mammals - 2,3,4 are too large

?max
