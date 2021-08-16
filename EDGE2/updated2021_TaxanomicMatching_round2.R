# matching Chondrichthyes data from fishbase with the trees downloaded from the shark website

# housekeeping

# Housekeeping
rm(list=ls())
graphics.off()
#setwd("C:\Users\oenon\AppData\Local\Packages\CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc\LocalState\rootfs\home\oslin")
getwd()
# set it as the linux machine - but have to use windows navigation as accessing through windows. GRIM!
setwd("../../../oenon/AppData/Local/Packages/CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc/LocalState/rootfs/home/oslin/EDGE2/")


load("../Data/ValidatedData.Rdata")

install.packages("ropensci/rfishbase")
update.packages("ropensci/rfishbase")

library(rfishbase)
library(phylobase)
library(phytools)
library(dplyr)

library(caper) #installed

library(data.table) #installed
library(geiger) #installed
library(pez)  #installed


###################

##### Import the three data sources we're going to be looking at

red <- read.csv( "Data/taxonomy.csv", stringsAsFactors = F)
#red <- red[-1]

tree <- read.csv("Data/PFC_taxonomy.csv", stringsAsFactors = F)
tree$scientificName <- tree$genus.species

fishbase <- as_tibble(rfishbase::load_taxa(
            server = "fishbase",
            version = "19.04",
            db = default_db(),
            collect = TRUE)) 

fishbase <- subset(fishbase, fishbase$Class == "Actinopterygii")

# nb - all binomial names have SPACES rather than "_". Leave like this for the comparison. 

# We are using the  fishbase species as the the ones that we trust
# we validate the tree species and the red list species

#########################################
#### Red & Fishbase 

fb.red  <- intersect(fishbase$Species, red$scientificName) # 17,999 species

# overlap between fishbase & species in the tree

fb.tree <- intersect(fishbase$Species, tree$scientificName) # 30,022 species

#don't validate - do one at a time
#fb.val <- validate_names(fishbase$Species) # somehow I have more species here? What!


###################################
# which species are in the tree but not in fishbase? 

Treeonly <- subset(tree, (tree$scientificName %in% fishbase$Species == F))
# 1494 species: check to see if any of them have been renamed - hopefully this reduces
# when we validate the tree species & check

Redonly <- subset(red, (red$scientificName %in% fishbase$Species == F))
# 1098 species are in the red list but not in FB - hopefully this reduces
# when we validate the tree species & check


############## validate names? 

chond.red$gs <- gsub("_", " ", chond.red$scientificName)
red.val <- validate_names(chond.red$gs) # of 1124 names there are 1112 valid names, i.e. 12 invalid ones!
red.val <- gsub(" ", "_", red.val)

# let's make it so it's obvious what is replaced for what
red.validated <- data.frame(red = unique(red$scientificName), validated = NA)
red.validated$red <- as.character(red.validated$red)

for (i in 1:length(red.validated$red)) {
  print(i)
  try(red.validated$validated[i] <- validate_names(red.validated$red[i]))
}


### now we look up the names that are so wrong that can't be automatically validated
# manual validations



# so now we can create a df that has the species names that align with fishbase, and get the red list assessment for them
# However, now we need to check what's going on with the tree species 

##########

# let's make it so it's obvious what is replaced for what

tree.validated <- data.frame(tree = unique(tree$scientificName), validated = NA)
tree.validated$tree <- as.character(tree.validated$tree)

for (i in 1:length(unique(tree.validated$tree))) {
  if (i %% 100 == 0) {print(i)}
  try(tree.validated$validated[i] <- validate_names(tree.validated$tree[i]))
}


### now let's manually check the names that aren't matching and fill them in 

save(red.validated, tree.validated, file = "../EDGE2/U2021_ValidatedData.Rdata")
save(red, tree, fishbase, file = "EDGEBaseData.Rdata")
save(Redonly, Treeonly, file = "EDGE_SppForMatching.Rdata")

#save(red.check, file = "RedCheck.Rdata")
#load("ValidatedData.Rdata")
#load("RedCheck.Rdata")

# how many NAs are there
length(red.validated$validated) - length(unique(red.validated$validated)) # 727

# how many NAs are there
length(tree.validated$validated) - length(unique(tree.validated$validated)) # 346 SAME!

length(unique(fishbase$Species)) - length(tree.validated$validated) # 1356 # also the same as before
# So there are <1.5k species that would need to be imputed



#  so, lets actually check 



######################### RESTART FROM HERE   

### gone through and fixed all of the species that I can. Let's now check the overlap

# starting with the red list 

red.check <- list()

for (i in 1:length(red.validated$red)) {
  if (is.na(red.validated$validated[i]) == T) {
    print(i)
    entry <- c(i, red.validated$red[i])
    red.check[[length(red.check)+ 1]] <- entry
  }
}

### now we have a list of all of the missing species and their index values, we investigate each species

# day 1 (9)
red.validated$validated[204] <- "Profundulus hildebrandi"
red.validated$validated[208] <- "Profundulus candalarius"
red.validated$validated[761] <- "Barbodes belinka"
red.validated$validated[790] <- "Acanthopsoides gracilentus"
red.validated$validated[883] <- "Acanthopsoides hapalias"
red.validated$validated[999] <- "Acanthopsoides gracilis"
red.validated$validated[1040] <- "Barbonymus schwanenfeldii"
red.validated$validated[1064] <- "Acanthopsoides delphax"
red.validated$validated[1154] <- "Helostoma temminckii"

# DAY 2 (c. 100)
red.validated$validated[1451] <- "Pachyurus squamipennis"
red.validated$validated[1475] <- "" # cant find
red.validated$validated[1508] <- "" # can't find
red.validated$validated[1542] <- "Parazacco fasciatus"
red.validated$validated[1560] <- "Laubuka siamensis"
red.validated$validated[1573] <- "Parazacco fasciatus"
red.validated$validated[1590] <- "Schistura incerta"
red.validated$validated[1616] <- "Pseudogobiopsis paludosus"
red.validated$validated[1618] <- "Macrognathus maculatus"
red.validated$validated[1837] <- "Bramocharax caballeroi"
red.validated$validated[1843] <- "Catostomus nebuliferus"
red.validated$validated[1968] <- "Profundulus labialis" # this one shares the same common name
red.validated$validated[1982] <- "Cichlasoma istlanum"
red.validated$validated[1985] <- "Cichlasoma trimaculatum"
red.validated$validated[2134] <- "Eviota rubra"
red.validated$validated[2495] <- "Gymnothorax dorsalis"
red.validated$validated[2932] <- "Poeciloconger kapala"
red.validated$validated[3118] <- "Cyprinion mhalensis"
red.validated$validated[3124] <- "Garra smarti"
red.validated$validated[] # 30-38 in list aren't described as species. Come back to them at the end
red.validated$validated[3847] <- "Padogobius nigricans"
red.validated$validated[4238] <- "Acanthocobitis botia"
red.validated$validated[4254] <- "Bangana ariza"
red.validated$validated[4324] <- "Clupisoma montana"
red.validated$validated[4336] <- "Pterocryptis indicus"
red.validated$validated[4364] <- "Garra rupecula"
red.validated$validated[4404] <- "Roeboides affinis"
red.validated$validated[4429] <- "Bryconamericus emperador"
red.validated$validated[4444] <- "Listrura camposi"
red.validated$validated[4493] <- "Trichomycterus stawiarski"
red.validated$validated[4521] <- "Steindachnerina corumbae"
red.validated$validated[4594] <- "Acanthocobitis rubidipinnis"
red.validated$validated[4662] <- "Chela khujairokensis"
red.validated$validated[4705] <- "Schistura acuticephalus"
red.validated$validated[4845] <- "Neoarius latirostris"
red.validated$validated[5014] <- "Akrokolioplax bicornis"
red.validated$validated[5231] <- "Laubuka fasciata"
red.validated$validated[5580] <- "Distichodus nefasch"
red.validated$validated[6734] <- # 57 - another unknown species
red.validated$validated[7663] <- "Suezichthys devisi"
red.validated$validated[] # 59 & 60 all unknown species
red.validated$validated[8048] <- "Cinetodus conorhynchus"
red.validated$validated[8095] # 62 also unknown species
red.validated$validated[8144] <- "Leptolebias opalescens"
red.validated$validated[8178] <- "Thalassophryne maculosa"
red.validated$validated[] # 65 - 86 are all unknown species
red.validated$validated[8212] <- "Laubuka caeruleostigmata"
red.validated$validated # 88-99 nebulous species
red.validated$validated[8625] <- "Leptolebias splendens"
red.validated$validated[8648] <- "Simpsonichthys constanciae"
red.validated$validated[] # 102 -110 unknown species 
red.validated$validated[10422] <- "Cobitis bilseli" # went with match for common name not old binomial
red.validated$validated # 112 - can't find. ANGRY pelican flounder
red.validated$validated [11151] <- # angry pelican flounder. Everywhere - but NOT fishbase!! 
red.validated$validated # 114 - another unknown species

# day three (c.40)
red.validated$validated[12058] <- "Crenidens indicus"
red.validated$validated[12080] <- "Aphyolebias peruensis"
red.validated$validated[12105] <- "Aphyolebias rubrocaudatus"
red.validated$validated[12112] <- "Aphyolebias obliquus"
red.validated$validated[12233] <- "Corydoras multimaculatus"
red.validated$validated[12275] <- "Orestias gymnota"
red.validated$validated[12289] <- "Aphyolebias wischmanni"
red.validated$validated[12579] <- "Maylandia nkhunguensis"
red.validated$validated[12689] <- "Aphyolebias claudiae"
red.validated$validated[12733] <- "Etroplus maculatus"
red.validated$validated[12959] <- #125 - can't find - but its DD
red.validated$validated[13013] <- "Chaetostoma aburrensis"
red.validated$validated[13625] <- "Pseudanthias fasciatus"
red.validated$validated[13846] <- "Monochirus atlanticus"
red.validated$validated[13892] <- "" #129 - can't find
red.validated$validated[13917] <- "" # 130 - also can't find
red.validated$validated[13923] <- "" # 131 can't find
red.validated$validated[14034] <- "Pomadasys crocro"
red.validated$validated[14146] <- "Telmatherina albolabiosus"
red.validated$validated[14164] <- "Oryzias bonneorum"
red.validated$validated[14197] <- "Barbodes platysoma"
red.validated$validated[14261] <- "Pseudogobiopsis festivus"
red.validated$validated[14283] <- NA # sometimes noted as synonym of  Lepidocephalus macrochir but these are listed sepeately in RL! Grr. 
red.validated$validated[14334] <- NA # 138 can't find URGH
red.validated$validated[14365] <- "Leiocassis aculeatus"
red.validated$validated[14366] <- "Leiocassis collinus"
red.validated$validated[14367] <- "Leiocassis tenebricus"
red.validated$validated[14402] <- "Pterocryptis taytayensis"
red.validated$validated[14435] <- "Pseudobagarius macronemus"
red.validated$validated[14485] <- "Betta breviobesus"
red.validated$validated[14588] <- NA # 145 questionably linked with Bathygobius fuscus, also known as Redigobius samberanoensis (not in fb)
red.validated$validated[14591] <- # 146 - undefined species
red.validated$validated[14646] <- # 147 - can't find 
red.validated$validated[14702] <- "" # 148 - can't find synonym
red.validated$validated[14706] # 149 - existings in eschmeyrs & RL not WoRMS or FB
red.validated$validated[14777] <- "Enteromius argenteus"
red.validated$validated[14780] <- "Pseudobarbus hospes"
red.validated$validated[14783] <- "Enteromius mattozi"
red.validated$validated[] # 153- not a defined species

# Day 4 (80)
red.validated$validated[14988] <- "" # confused with 'shorthead anchovy' - but now with seperate spp acknowledgement
red.validated$validated[14989] <- "" # also not able to be found
red.validated$validated[14991] <- "" #  156 - can't find (also DD)
red.validated$validated[14992] <- "" # 157 - can't find, DD
red.validated$validated[15000] <- "" # 158 - can't find, DD
red.validated$validated[15002] <- "" # 159 - can't find, DD
red.validated$validated[15005] <- "" # 160 FB no match, Esch says Macrocephenchelys brevirostris but seperaetely listed
red.validated$validated[15007] <- "" # 161 actually Bathycongrus macroporis, but still no match
red.validated$validated[] # 162 - 173 not specific species
red.validated$validated[15051] <- "" # 174 can't find - DD 
red.validated$validated[15067] <- "Bathyclupea megaceps"
red.validated$validated[15068] <- "Bathyclupea argentea"
red.validated$validated[15083] <- "" # 177 said to be synonymous with Etrumeus teres however found in different places
red.validated$validated[15128] <- "" # 178 no match
red.validated$validated[15160] <- "Oligolepis keiensis"
red.validated$validated[15945] <- "Pseudocaranx dentex"
red.validated$validated[15952] <- "" # 181 - can't find, DD
red.validated$validated[15963] <- "Hemigrammocypris rasborella"
red.validated$validated[15968] <- "Pseudobagrus aurantiacus"
red.validated$validated[15969] <- "" # 184 - NA, assocaited with "Pseudobagrus aurantiacus" howver that is allocated to Tachysurus tokiensis
red.validated$validated[16043] <- "Coryogalops sordida"
red.validated$validated[16062] <- "Laciris pelagicus"
red.validated$validated[16080] <- "" # 187 - no synonyms found 
red.validated$validated[16081] <- "" # 188 - no synonyms found
red.validated$validated[16099] <- "" # 189 - no synonyms found, DD
red.validated$validated[16141] <- "Micracanthicus vandragti"
red.validated$validated[16151] <- "" # 191 - no synonyms found
red.validated$validated[16170] <- "Oxyurichthys microlepis" # already in system 
red.validated$validated[16171] <- "Oxyurichthys petersii" # also already accounted for
red.validated$validated[16183] <- "Trimma filamentosus" # unsure why didn't get picked up
red.validated$validated[16199] <- "" # 195 - no synonym, no match , DD
red.validated$validated[16200] <- "Epinephelus chlorostigma" # already in db
red.validated$validated[16214] <- "" # no synonym
red.validated$validated[16216] <- "Taeniolethrinops praeorbitalis" # already in db
red.validated$validated[16219] <- "Maylandia lanisticola"
red.validated$validated[16268] <- "" # can't find - DD -  taxanomic review of this whole fam Family Ateleopodidae needed
red.validated$validated[16272] <- "Ijimaia fowleri"
red.validated$validated[16296] <- "" # 202 - no synonym
red.validated$validated[16302] <- "" # 203 - no synonym
red.validated$validated[16323] <- "" # 204 - no synonym
red.validated$validated[16346] <- "Lestrolepis luetkeni"
red.validated$validated[16355] <- "Notolepis coatsi"
red.validated$validated[16385] <- "Benthalbella macropinna"
red.validated$validated[16395] <- "" # 208 - no synonym & DD
red.validated$validated[] # 209 - no synonym & DD
red.validated$validated[16407] <- "" # 210 - no synonym 
red.validated$validated[16418] <- "" # 211 - no synonym
red.validated$validated[16485] <- "" # 212 - suggested Lophotus lacepede, but already in list
red.validated$validated[16565] <- "" # 213 - no synonym, unfortunately
red.validated$validated[16566] <- "" # 214 - no synonym, split from Galaxias olidus
red.validated$validated[16568] <- "" # 215 - no synonym, split from Galaxias olidus
red.validated$validated[16569] <- "" # 216 - no synonym, split from Galaxias olidus
red.validated$validated[16570] <- "" # 217 - no synonym, split from Galaxias olidus
red.validated$validated[16571] <- "" # 218 - no synonym, split from Galaxias olidus
red.validated$validated[16572] <- "" # 219 - no synonym, split from Galaxias olidus
red.validated$validated[16573] <- "" # 220 - no synonym, split from Galaxias olidus
red.validated$validated[16575] <- "" # 221 - no synonym, split from Galaxias olidus
red.validated$validated[16577] <- "" # 222 - no synonym, split from Galaxias olidus
red.validated$validated[16578] <- "" # 223 - no synonym, split from Galaxias olidus
red.validated$validated[16579] <- "" # 224 - no synonym, split from Galaxias olidus
red.validated$validated[] # 225 - 229 ill defined species
red.validated$validated[16689] <- "" # 230 - no synonym
red.validated$validated[16704] <- "Parachaeturichthys ocellatus"
red.validated$validated[16741] <- "" # 232 - no synonm, DD
red.validated$validated[16781] <- "Engyprosopon sechellensis"

# Day five (c.80)

red.validated$validated[16984] <- "Enteromius quadrilineatus"
red.validated$validated[17018] <- "Crossocheilus latius"
red.validated$validated[17021] <- "Chelon planiceps"
red.validated$validated[] # 237 - 241 - undefined species
red.validated$validated[17043] <- "" # 242 - no synonym, new species.Used to be associated with Cairnsichthys rhombosomoides, however distinct population
red.validated$validated[17045] <- "Guyu wujalwujalensis"
red.validated$validated[17072] <- "Macquaria novemaculeata"
red.validated$validated[17082] <- "" # 245 - DD - no synonym & validity questioned
red.validated$validated[17107] <- "" # 256 - DD - ambiguous synonym of Lampanyctus festivus but listed seperately
red.validated$validated[17108] <- "Nannobrachium bristori" 
red.validated$validated[17109] <- "Nannobrachium crypticum"
red.validated$validated[17111] <- "Lampanyctus crocodilus" # already in db - invalid name before
red.validated$validated[17112] <- "Nannobrachium hawaiiensis"
red.validated$validated[17115] <- "Nannobrachium indicum"
red.validated$validated[17121] <- "Nannobrachium phyllisae"
red.validated$validated[17129] <- "Triphoturus nigrescens" # already in db
red.validated$validated[17208] <- "Catostomus plebeius"
red.validated$validated[17248] <- "Myctophum brachygnathum"
red.validated$validated[17273] <- "" # 256 - DD, no synonym
red.validated$validated[17318] <- "" # 257 - used to be Nothobranchius rosenstocki however recently seperated
red.validated$validated[17319] <- "" # 258 - no synonym
red.validated$validated[17321] <- "" # 259 - no synonyn, recent split from Nothobranchius lucius & melanospilus
red.validated$validated[17346] <- "Astyanax fasciatus"
red.validated$validated[17347] <- "" # this taxonomy is highly confused - so many synonyms, all of different species..
red.validated$validated[17349] <- "" # ditto above
red.validated$validated[17367] <- "Petroleuciscus kurui" # DD though
red.validated$validated[17369] <- "Profundulus mixtlanesis"
red.validated$validated[17371] <- "" # DD
red.validated$validated[17385] <- "" # DD 
red.validated$validated[17386] <- "" # synonym is already in db
red.validated$validated[17392] <- "Heterotilapia buttikoferi"
red.validated$validated[17400] <- "Coptodon rheophila"
red.validated$validated[17427] <- "" # no synonym, DD
red.validated$validated[17428] <- "" # no synonym, DD
red.validated$validated[17429] <- "" # synonym already in db, DD anyway
red.validated$validated[17432] <- "" # no clear synonym, DD
red.validated$validated[17439] <- "" # no synonym
red.validated$validated[17467] <- "" # no clear synonym, DD
red.validated$validated[17468] <- "" # no synonym
red.validated$validated[17469] <- "" # neither name nor synonym found 
red.validated$validated[17470] <- "" # no synonym
red.validated$validated[17471] <- "" # no synonym, DD
red.validated$validated[17472] <- "" # no synonym, DD
red.validated$validated[17473] <- "" # no synonym,
red.validated$validated[17474] <- "" # no synonym, DD
red.validated$validated[17475] <- "" # no synonym, DD
red.validated$validated[17476] <- "" # no synonym, DD
red.validated$validated[17477] <- "" # no synonym, DD
red.validated$validated[17478] <- "" # no synonym, DD
red.validated$validated[17480] <- "Lipogramma haberi"
red.validated$validated[17478] <- "" #  no synonym, DD  
red.validated$validated[17505] <- "Melamphaes longivelis" # already in db, dd 
red.validated$validated[17506] <- ""  # synonym alrady in db, this dd other rated
red.validated$validated[17526] <- "" # Esch synonym already in db:  Astronesthes niger
red.validated$validated[17545] <- ""# no synonym, dd
red.validated$validated[17550] <- "" # no synonym, dd
red.validated$validated[17560] <- "" # no synonym, dd
red.validated$validated[17566] <- "" # no synonym, dd
red.validated$validated[17626] <- "Barbus lorteti" 
red.validated$validated[17634] <- "Labeobarbus rhinoceros" # same name - false flag
red.validated$validated[17642] <- "" # synonym already in db, Luciobarbus callensis
red.validated$validated[17643] <- "Pseudophoxinus callensis"
red.validated$validated[17646] <- "" # here synonym already in db. Luciobarbus callensis
red.validated$validated[17649] <- "" # again synonym already in db. Luciobarbus callensis
red.validated$validated[17661] <- "Squalius aphipsi"
red.validated$validated[17663] <- "Pomatoschistus canestrinii"
red.validated$validated[17664] <- "Pomatoschistus montenegrensis"
red.validated$validated[17665] <- "Knipowitschia punctatissima"
red.validated$validated[17681] <- "" # synonym already in db
red.validated$validated[17687] <- "" # no valid synonym
red.validated$validated[17715] <- "Fiordichthys paxtoni"
red.validated$validated[17716] <- "Fiordichthys slartibartfasti"
red.validated$validated[17725] <- "Diplacanthopoma brunnea"

# day 6
red.validated$validated[17742] <- "Cataetyx niki"
red.validated$validated[17743] <- "Parasciadonus pauciradiatus"
red.validated$validated[17751] <- "Aphyonus bolini"
red.validated$validated[17752] <- "Aphyonus brevidorsalis"
red.validated$validated[17753] <- "" # 315 no synonym, DD
red.validated$validated[17754] <- "Barathronus solomonensis"
red.validated$validated[17764] <- "" # 317, no synonym, DD
red.validated$validated[17798] <- "" # 318, no synonym, DD
red.validated$validated[17809] <- "Apagesoma australis"
red.validated$validated[17810] <- "" # 320, no synonym, DD
red.validated$validated[17811] <- "" # 321, no synonym, DD
red.validated$validated[17846] <- "" # 322, no synonym, DD
red.validated$validated[17852] <- "" # 323, no synonym
red.validated$validated[17856] <- "" # 324, no synonym
red.validated$validated[17862] <- "" # 325. no synonym, DD  
red.validated$validated[17872] <- "" # 326, no synonym
red.validated$validated[17923] <- "" # 327, no synonym
red.validated$validated[17946] <- "" # 328. no synonym, DD
red.validated$validated[17953] <- "" # 329, no synonym
red.validated$validated[17956] <- "" # 330, no synonym, DD
red.validated$validated[18024] <- "" # 331, no synonym
red.validated$validated[18033] <- "Melamphaes papavereus" # found the sp Melamphaes papavereus in fishbase, but not coming up in db
red.validated$validated[18036] <- "" # 333, no synonym
red.validated$validated[18037] <- "Melamphaes succedaneus"
red.validated$validated[18065] <- "" # no synonym, DD
red.validated$validated[18071] <- "" # no synonym
red.validated$validated[18072] <- "" # no synonym
red.validated$validated[18074] <- "Rhamdia guatemalensis"
red.validated$validated[18075] <- "" # 339 - no synonym
red.validated$validated[18076] <- "Glyptothorax platypogon" # found in fb, but saying 'not in download'. Odd!
red.validated$validated[18077] <- "" # 341, no synonym
red.validated$validated[18082] <- "Chaunax abei" # already in db
red.validated$validated[18087] <- "" # 343, no synonym
red.validated$validated[18092] <- "" # 344, no synonym
red.validated$validated[18118] <- "" # 345, no synonym, DD
red.validated$validated[18119] <- "" # 346. no synonym
red.validated$validated[18126] <- "" # 347, no synonym, DD
red.validated$validated[18189] <- "" # 348, no synonym, DD
red.validated$validated[18192] <- "" # 349, no synonym, DD
red.validated$validated[18203] <- "" # 350, no synonym, DD
red.validated$validated[18207] <- "" # 351, no synonym, DD
red.validated$validated[18210] <- "" # 352, no synonym, DD
red.validated$validated[18212] <- "" # 353, no synonym, DD
red.validated$validated[18216] <- "" # 354. no synonym, DD
red.validated$validated[18256] <- "Hypoplectrus liberte" # claims on in FB, but IN FB! So will use. 
red.validated$validated[18308] <- "" # 357, no sysonym, DD
red.validated$validated[18330] <- "" # 358, no synonym, DD
red.validated$validated[18333] <- "" # 359. no synonym, DD
red.validated$validated[18338] <- "Stomias longibarbatus" # already in db
red.validated$validated[18349] <- "Starksia splendens" # claims not in FB, but there!
red.validated$validated[18377] <- "" # 362, no synonym, DD
red.validated$validated[18380] <- "" # 363, no synonym, DD
red.validated$validated[18425] <- "Spottobrotula amaculata" 
red.validated$validated[18429] <- "" # 365, no synonym, DD
red.validated$validated[18431] <- "" # 366, no synonym, DD
red.validated$validated[18432] <- "" # 367, no synonym, DD
red.validated$validated[18441] <- "" # 368, no synonym, DD
red.validated$validated[18442] <- "" # 369, no synonym, DD
red.validated$validated[18443] <- "Enneanectes flavus" # not in FB but in fb!
red.validated$validated[18450] <- "" # 371, no synonym, DD
red.validated$validated[18475] <- "" # 372, no synonym, DD
red.validated$validated[18465] <- "" # 373, no synonym, DD
red.validated$validated[18517] <- "Sundasalanx malleti"
red.validated$validated[18555] <- "" # 375. no synonym, DD
red.validated$validated[18558] <- "" # 376, no synonm
red.validated$validated[18565] <- "" # Psettina profunda already in red list 
red.validated$validated[18582] <- "" # 378, no synonym, DD
red.validated$validated[18587] <- "Engyprosopon annulatus"
red.validated$validated[18589] <- "" # 380, no synonym, DD
red.validated$validated[18593] <- "Engyprosopon kushimotoensis"
red.validated$validated[18598] <- "" # 382, no synonym, DD
red.validated$validated[18599] <- "" # Arnoglossus aspilos already in db 
red.validated$validated[18627] <- "" # Gobiodon rivulatus already in db
red.validated$validated[18630] <- "" # already in db
red.validated$validated[] # 386 - 396 not defined species
red.validated$validated[18645] <- "" # Cynoglossus dollfusi already in db
red.validated$validated[18652] <- ""  # 398, no synonym, DD
red.validated$validated[18655] <- "" # 399, no synonym
red.validated$validated[18687] <- "" # 400 neither name nor synonym recognised
red.validated$validated[18773] <- "Maylandia cyneusmarginata"
red.validated$validated[18794] <- "Trematocranus pachychilus" # says no in FB but in FB!
red.validated$validated[18814] <- "Eviota lentiginosa" # says not in FB but in FB!
red.validated$validated[18857] <- ""  # 404, no synonym, 
red.validated$validated[18861] <- "Acanthocobitis urophthalmus"
red.validated$validated[18863] <- "" # Schistura sertata already in db
red.validated$validated[18864] <- ""# no synonym
red.validated$validated[18884] <- "" # 408 - 418 not exact species
red.validated$validated[18905] <- "Laubuka insularis" # name recognised
red.validated$validated[18906] <- "Laubuka lankensis" # name recognised
red.validated$validated[18907] <- "Laubuka ruhuna" # name recognised
red.validated$validated[18908] <- "Laubuka varuna"  # name recognised
red.validated$validated[19022] <- "Protomelas spilopterus"
red.validated$validated[19037] <- "" # 424, no synonym
red.validated$validated[19038] <- "" # 425, no synonym
red.validated$validated[19039] <- "" # 426, no synonym
red.validated$validated[19404] <- "" # 427, no synonym
red.validated$validated[19041] <- "" # 428, no synonym
red.validated$validated[19042] <- "" # 429, no synonym
red.validated$validated[19051] <- "" # 430, no synonym
red.validated$validated[19058] <- "" # 431, no synonym
red.validated$validated[19061] <- "" # 432, no synonym
red.validated$validated[19062] <- "" # 433, no synonym
red.validated$validated[19063] <- "" # 434, no synonym
red.validated$validated[19064] <- "Sardinella pacifica" # on fb bur coming up false
red.validated$validated[19065] <- "" # 436, no synonym
red.validated$validated[19066] <- "" # 437, no synonym
red.validated$validated[19067] <- "" # 438 synonym already in db Bairdiella ronchus
red.validated$validated[19068] <- "" # 439, no synonym
red.validated$validated[19069] <- "" # 440, no synonym
red.validated$validated[19075] <- "" # 441, no synonym
red.validated$validated[19077] <- "" # 442, no synonym
red.validated$validated[19078] <- "" # 443, synonym already in db Enteromius apleurogramma
red.validated$validated[19081] <- "Sciaena wieneri"
red.validated$validated[19090] <- "" # 445, synonym already in db Kentrocapros aculeatus


red.validated$validated[x] %in% fishbase$Species # need to return true
# also worth checking whether it's already in the red list? 
red.validated$validated[x] %in% red$scientificName # ideally return false

y <- 445
red.check[[y]]
x <- as.integer(red.check[[y]][1])

# checker
red.validated$validated[x] %in% fishbase$Species

save(red.validated, tree.validated, total.add, file = "UpdatedValidSpp.Rdata")
# next one

length(red.check)

###
add.to.fb <- c("Glyptothorax platypogon", "Hypoplectrus liberte", "Starksia splendens", 
               "Enneanectes flavus", "Trematocranus pachychilus", "Eviota lentiginosa", 
               "Sardinella pacifica")

#####

##### now we do the same with tree.validated 

tree.check <- list()

for (i in 1:length(tree.validated$tree)) {
  if (is.na(tree.validated$validated[i]) == T) {
    print(i)
    entry <- c(i, tree.validated$tree[i])
    tree.check[[length(tree.check)+ 1]] <- entry
  }
}

length(tree.check) # only 66 to go! 

red.na <- subset(red.validated, is.na(red.validated$validated) == T)
tree.na <- subset(tree.validated, is.na(tree.validated$validated) == T)

na.overlap <- intersect(red.na$red, tree.na$tree)
# now claiming to be zero. how weird. I guess we'll start from scratch with the tree.check then. 


## checks
tree.validated$validated[589] <- "Leptocephalus bellottii"
tree.validated$validated[1030] <- "Ijimaia fowleri"
tree.validated$validated[1773] <- "Badis triocellus"
tree.validated$validated[3270] <- "" # synonym already in db
tree.validated$validated[3339] <- "Kyphosus azureus"
tree.validated$validated[5064] <- "Sagittalarva inornatus"
tree.validated$validated[5156] <- "Suezichthys devisi"
tree.validated$validated[5211] <- "Xyrichtys javanicus"
tree.validated$validated[5744] <- "" # referring to a third Enigmapercis spp not in fb
tree.validated$validated[6642] <- "" # no synonym
tree.validated$validated[7338] <- "" # no synonym
tree.validated$validated[7414] <- "" # synonym already in db
tree.validated$validated[7484] <- "" # no synonym
tree.validated$validated[7485] <- "Nothonotus microlepidus"
tree.validated$validated[7536] <- "Nothonotus starnesi"
tree.validated$validated[7583] <- "" #  16, not a described species
tree.validated$validated[8329] <- "Liopropoma danae"
tree.validated$validated[9213] <- "" # 18 no synonym
tree.validated$validated[9361] <- "" # 19 no synonym
tree.validated$validated[9571] <- "Astroscopus y-graecum"
tree.validated$validated[10571] <- "Gobiopterus smithi"
tree.validated$validated[10558] <- "Gobius tropicus"
tree.validated$validated[11118] <- "Redigobius balteatops"
tree.validated$validated[11396] <- "Tigrigobius digueti"
tree.validated$validated[11667] <- "Apogon latus"
tree.validated$validated[13405] <- "" # 26, "Chirostoma estor" synonym already in db, 
tree.validated$validated[14432] <- "" # 27, "Simpsonichthys bokermanni" already in db
tree.validated$validated[14433] <- "" # 28, "Simpsonichthys constanciae", already in db
tree.validated$validated[14795] <- "Lacustricola bukobanus"
tree.validated$validated[16393] <- "" # not defined Astronotus species 
tree.validated$validated[16620] <- "Coptodon bythobates"
tree.validated$validated[16621] <- "Coptodon camerunensis"
tree.validated$validated[17082] <- "" # not a defined species - there are many Hemichromis sp 
tree.validated$validated[17299] <- "Trichromis salvini"
tree.validated$validated[17585] <- "" # not a defined species
tree.validated$validated[18414] <- "" # no matching species
tree.validated$validated[18881] <- "" # synonym already in db
tree.validated$validated[20776] <- "" # 39 - not a defined species
tree.validated$validated[21262] <- "Brycinus batesii"
tree.validated$validated[21290] <- "Brycinus schoutedeni"
tree.validated$validated[21334] <- "" # synonym already in db (Alestopetersius compressus)
tree.validated$validated[21339] <- "" # synonym already in db (Alestopetersius nigropterus)
tree.validated$validated[21657] <- "" # synonym already in db (Astyanax jordani)
tree.validated$validated[21721] <- "" # not a defined species
tree.validated$validated[22538] <- "Roeboides affinis"
tree.validated$validated[23296] <- "" # not a defined species
tree.validated$validated[23648] <- "" # no synonym, not Corydoras multimaculatus as in tehre seperately
tree.validated$validated[23727] <- "" # synonym (Corydoras kronei) not in fb
tree.validated$validated[25113] <- "" # synonym (Netuma proxima) already in list
tree.validated$validated[25219] <- "" # synonmy (Ariopsis coatesi) not recognised
tree.validated$validated[25220] <- "" # synnym (Neoarius latirostris) already in fb
tree.validated$validated[25222] <- "" # synonym (Neoarius velutinus) already in fb
tree.validated$validated[26961] <- "" # synonym (Glyptothorax zanaensis) already in fb
tree.validated$validated[27341] <- "" # synonym (Yasuhikotakia lecontei) already in fb
tree.validated$validated[27523] <- "" # not a defined species
tree.validated$validated[27544] <- "" # synonym already in fb (Cobitis pacifica)
tree.validated$validated[27548] <- "" # synonym a) already in fb (Cobitis multifasciata) synonym b) Cobitis multifasciata isn't used
tree.validated$validated[28965] <- "Bangana behri"
tree.validated$validated[29515] <- "" # synon ym alraedy validated list 
tree.validated$validated[29851] <- "" # synonym (Pseudophoxinus callensis) already in validated list
tree.validated$validated[30470] <- "Systomus orphoides"
tree.validated$validated[30706] <- "" # synonym (Lepidocephalichthys guntea) in validated list
tree.validated$validated[31495] <- "" # synonym (Scleropages formosus) in validated list
tree.validated$validated[31499] <- "" # ambig synonym of S. formosus also - not recognised alone
tree.validated$validated[31501] <- "" # also ambig synonym


tree.validated$validated[x] %in% fishbase$Species # need to return true
# also worth checking whether it's already in the red list? 
tree.validated$validated[x] %in% tree$scientificName # ideally return false


y <- 67
tree.check[[y]]
x <- as.integer(tree.check[[y]][1])

# checker
tree.validated$validated[x] %in% fishbase$Species


tree.to.fb <-c("Gobiopterus smithi")
##### pick up with matching the tree spp tomorrow, then create a new master list 
# (inc. the few from fb I need to add in) and do the imputations - check with Rikki
# if there is a source of 1000 Actinopt trees, or whether I just use the 100 I have randomly. 

total.add <- c(tree.to.fb, add.to.fb)



############################# now create the master species list 

fml <- data.frame(Species = c(fishbase$Species, total.add), 
                               RLScore = NA, TreeSpp = NA, Impute = NA)

fml$Species <- as.character(fml$Species)

to.impute <- subset(fml, (fml$Species %in% tree.validated$validated == F)) # 1691 species

red$validated <- red.validated$validated

for (i in 1:length(fml$Species)) {
  if (i %% 100 == 0) {
    print(i)
  }
  ref <- match(fml$Species[i], red$validated)
  if (length(ref) != 0) {
    fml$RLScore[i] <- red$redlistCategory[ref]
  }
}

###

for (i in 1:length(fml$Species)) {
  if (i %% 100 == 0) {
    print(i)
  }
  ref <- match(fml$Species[i], tree.validated$validated)
  if (length(ref) != 0) {
    fml$TreeSpp[i] <- tree.validated$tree[ref]
  }
}

###

fml$Impute <- 1

for (i in 1:length(fml$Species)) {
  if (fml$Species[i] %in% tree.validated$validated) {
    fml$Impute[i] <- 0
  }
}

sum(fml$Impute) # 1691 species to impute. Phew! 


##############

# facts & figures - comparing the overlaps 

length(intersect(tree.validated$validated, fishbase$Species)) #31187
length(tree.validated$validated) # 31516
length(unique(tree.validated$validated)) # 31190 - worth noting that we still have 10 entries that are NA

####

save(fml, file = "FishAllSpecies.Rdata")

