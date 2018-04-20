library(CBDD)
source("~/GITCLONES/cbdd/load.R")
setwd("~/GITCLONES/2017-12-07-kl2_kl4_upstream_analysis/system_biology")


########################################################################################################
#Step 0:
#You could load these data and use it directly instead of run next steps
#
#
load("jp_work.RData")
load("ReactomePathwayList.Rdata")
load("klf2_klf4_upstream_testing.RData")
load("FIs_network.Rdata")
########################################################################################################




########################################################################################################
#Step 1:
################################################################
#All the files use here were downloaded from reactome
#https://reactome.org/download-data/
#Here I used Functional Interactions derived from Reactome, and other pathway and interaction databases(version 2016)
#and PE Identifier mapping files(NCBI to All pathways)
################################################################

library(data.table)
a <- fread("FIsInGene_022717_with_annotations.txt",sep = "\t",header = T) 
#filter wrong direction
Step1 <- a[!startsWith(a$Direction,"-"),]
#filter double direction
Step2 <- Step1[grep("[<][-][>]",Step1$Direction),]
Step3 <- Step1[grep("[<][-][|]",Step1$Direction),]
Step4 <- Step1[grep("[|][-][>]",Step1$Direction),]
Step5 <- Step1[grep("[|][-][|]",Step1$Direction),]
Step2$Direction <- "->"
Step2_1 <- Step2
colnames(Step2_1) <- c("Gene2",      "Gene1"  ,    "Annotation", "Direction"  ,"Score"   )
Step5$Direction <- "-|"
Step5_1 <- Step5
colnames(Step5_1) <- c("Gene2",      "Gene1"  ,    "Annotation", "Direction"  ,"Score"   ) 
Step3$Direction <- "-|"
Step3_1 <- Step3
colnames(Step3_1) <- c("Gene2",      "Gene1"  ,    "Annotation", "Direction"  ,"Score"   ) 
Step3_1$Direction <- "->"
Step4$Direction <- "->"
Step4_1 <- Step4
colnames(Step4_1) <- c("Gene2",      "Gene1"  ,    "Annotation", "Direction"  ,"Score"   ) 
Step4_1$Direction <- "-|"

#first reverse non "<->"
Step6 <- Step1[-c(grep("[<][-][>]",Step1$Direction),grep("[<][-][|]",Step1$Direction),grep("[|][-][>]",Step1$Direction),grep("[|][-][|]",Step1$Direction)),]
colnames(Step6) <- c("Gene2",      "Gene1"  ,    "Annotation", "Direction"  ,"Score"   )
#replace"<-"," |-"
Step6$Direction <- gsub("<-","->",Step6$Direction)
Step6$Direction <- gsub("[|][-]","-|",Step6$Direction)

################collect and combine
wellformat <- a[startsWith(a$Direction,"-"),]
wellformat$mechanism <- 0
wellformat$mechanism[grep("expression regula",wellformat$Annotation)] <- 1

draft <- rbind(Step2,Step2_1,Step3,Step3_1,Step4,Step4_1,Step5,Step5_1,Step6)
draft$mechanism <- 1
draft <- rbind(wellformat,draft)
draft$Direction <- gsub("[-][>]","1",draft$Direction )
draft$Direction <- gsub("[-][|]","-1",draft$Direction )
draft$Direction <- gsub("[-]$","0",draft$Direction )
colnames(draft) <- c("Gene1",      "Gene2"  ,    "mechanism_name", "edgetype"  ,"Trust" ,"mechanism"  )
draft$directed <- "0"
draft$directed[draft$edgetype != 0] <- 1
draft$edgetype <- as.numeric(draft$edgetype)
homoFIs <- draft
save(homoFIs,file = "FIs_network.Rdata")

###build pathways
library(data.table)
library(stringr)
a <- fread("NCBI2Reactome_PE_All_Levels.txt",sep = "\t",header = F) 
a1 <- data.frame(str_split_fixed(a$V3, "[ ]", 2))
a2 <- a
a2$V3 <- gsub("^[ ]","",a2$V3)
a2 <- a2[,-c(5,7)]
a2$localization <- a1$X2
colnames(a2) <- c("NCBI","RPESI","RPE","Pathway","PathwayName","speices","localization")
a2$RPE <- a1$X1
ReacPathway <- a2 
#remove parenthesis
ReacPathway$RPE <- gsub("\\([^\\)]+\\)","",ReacPathway$RPE)
#filter by homo sapiens
ReacPathway <- ReacPathway[grep("Homo sapiens",ReacPathway$speices),]
#remove -
non_nu <- ReacPathway[is.na(as.numeric(gsub(".*\\-","",ReacPathway$RPE))),]
non_nu$RPE <- factor(non_nu$RPE)
non_nu$RPE <- gsub(".*\\-","",ReacPathway$RPE)
#View(data.frame(non_nu$RPE)) check for comma
non_nu <- non_nu[-grep("\\,",non_nu$RPE),]
nu <- ReacPathway[!is.na(as.numeric(gsub(".*\\-","",ReacPathway$RPE))),]
nu$RPE <- gsub('.*\\,', '',nu$RPE)
nu$RPE <- gsub(".*\\-(.*)\\-.*", "\\1", nu$RPE)
nu$RPE <-  gsub("\\-.*","",nu$RPE)

full <- rbind(nu,non_nu)

# oo <- data.frame( name="names")
# for(i in 1:dim(ReacPathway)[1]){
# if(is.na(as.numeric(gsub(".*\\-","",ReacPathway$RPE[i])gsub(".*\\-","",ReacPathway$RPE[i])gsub(".*\\-","",ReacPathway$RPE[i]))) == T){
# #ReacPathway$RPE[i] <- 
#   tt <- data.frame(gsub(".*\\-","",ReacPathway$RPE[i]))
#                   colnames(tt) <- "name"
#   oo <- rbind(oo,tt)
# }else{
#   tt <- data.frame(ReacPathway$RPE[i])
#   colnames(tt) <- "name"
#   oo <- rbind(oo,tt)
# }
# }
ReacPathwayList <- split( full, f = full$PathwayName)


#write a loop to seperate the list to data frame and clean out cut-offs (still log tpm)
name_list <- names(ReacPathwayList)
int_list <- list()
for(i in 1:length(name_list)){
  #a<- gsub("\\s","",name_list[i])
  a <-list(as.character(ReacPathwayList[[i]]$RPE))
  names(a) <- name_list[i]
  int_list <- append(a,int_list)
  
}

temp_list <- list()
temp_empty <- list()
saved_int <- int_list
for( i in 1:length(int_list)){
  int_list[[i]] <- na.omit(unname(convert[int_list[[i]]]))
  
}
EMhomo <- homoFIs
EMhomo$Gene1 <- unname(convert[EMhomo$Gene1])
EMhomo$Gene2 <- unname(convert[EMhomo$Gene2])
EMhomo <- na.omit(EMhomo)
for (i in 1:length(int_list)){
  
  temp <- EMhomo[(as.character(EMhomo$Gene1) %in% as.character(unlist(unname(int_list[i]))) & as.character(EMhomo$Gene2) %in% as.character(unlist(unname(int_list[i])))),]
  if(dim(temp)[1] != 0){
    
    temp$link_id <- rep(0,dim(temp)[1])
    
    temp <- list(temp)
    names(temp) <-names(int_list[i])
    temp_list <- append(temp,temp_list)
  }else{
    print(paste(names(int_list[i]),"is empty"))
    temp_empty <- append(temp_empty,list(names(int_list[i])))
  }
  
  
}

pathway_list <- lapply(temp_list,as.data.frame)
save(pathway_list,file = "~/Desktop/ReactomePathwayList.Rdata")



########################################################################################################






########################################################################################################
#Step 3: test KLF2 & KLF4 using these two networks.
# NOTE: Since we do not have mechanism infor for reactome networks, we set it default to 1
#
#################################################################

mini_network <- network
mini_pathway <- pathways
reac_network <- as.data.frame(homoFIs)
reac_pathway <- pathway_list

mini_to_reac1 <- as.character(minibase.id.converter$HGNC_SYMBOL)
names(mini_to_reac1) <- as.character(minibase.id.converter$NETWORK_OBJECT_ID)
#mini_network <- na.omit(mini_network)
#mini_network$network_object_id1 <- unname(mini_to_reac1[mini_network$network_object_id1])
#mini_network$network_object_id2 <- unname(mini_to_reac1[mini_network$network_object_id2])
persudo_logfc <- data.frame( gene <- c("KLF2","KLF4"),
                             logfc <- c(10,10)
)
mini_logfc <- data.frame(ID <- c("-2007163176","-2066067575"),
                         logfc <- c(10,10)
)
try <- reac_network
reac_network <- try[,c(1,2,4,6,3,5,7)]
colnames(reac_network) <-c("network_object_id1", "network_object_id2" ,"edgetype"    ,       "mechanism"  ,       
                          "mechanism_name",     "trust","directed" )
reac_network$mechanism <-1
test_mini <- SigNet(mini_logfc,mini_network)
test_reac <- SigNet(persudo_logfc,reac_network)
test_mini$gene <- unname(mini_to_reac1[test_mini$hypothesis])
#############try other functions
#2. causal reasoning 
# test2_mini <- causalReasoning(mini_logfc,mini_network,scoreThreshold = 0, enrichment = F)
# test2_reac <- causalReasoning(persudo_logfc,reac_network, scoreThreshold = 0,enrichment = F)
# test2_mini$hypothesis <- unname(mini_to_reac1[test2_mini$hypothesis])
# #3. neighborhood scoring
# test3_mini <- neighborhoodScoring(mini_logfc,mini_network)
# test3_reac <- neighborhoodScoring(persudo_logfc,reac_network)
# 
# #4. guilt by asso
# test4_mini <- guiltByAssociation(mini_logfc,mini_network)
# test4_reac <- guiltByAssociation(persudo_logfc,reac_network)
# 
# #5. interconnec
# 
# test5_mini <- interconnectivity(mini_logfc,mini_network)
# test5_reac <- interconnectivity(persudo_logfc,reac_network)

#6. network propagation
test6_mini <- networkPropagation(mini_logfc,mini_network,directed = T)
test6_reac <- networkPropagation(persudo_logfc,reac_network,directed = T)

reverse_mininetwork <- mini_network
reverse_mininetwork$edgetype<- gsub("-1","9",reverse_mininetwork$edgetype)
reverse_mininetwork$edgetype<- gsub("1","-1",reverse_mininetwork$edgetype)
reverse_mininetwork$edgetype<- gsub("9","1",reverse_mininetwork$edgetype)


test6_mini_rev <- networkPropagation(mini_logfc,reverse_mininetwork,directed = T)
test6_reac_rev <- networkPropagation(persudo_logfc,reac_network,directed = T)
# 
# #7. randomwalk
# test7_mini <- randomWalk(mini_logfc,mini_network)
# test7_reac <- randomWalk(persudo_logfc,reac_network)
# #rev
# test7_mini_rev <- networkPropagation(mini_logfc,reverse_mininetwork,directed = T)
# 
# 
# 
# #8. topnet
# test8_km_mini$gene <- unname(mini_to_reac1[test8_km_mini$Node_identifier])
# 
# levels <- numeric()
# n <- dim(test8_km_mini)[1]
# levels <- rep(0,as.numeric(n))
# for (i in 1:length(levels)){
#   a <-test8_km_mini$gene[i]
#   for(j in 1:10){
#     b <- up_list[[j]]
#     if( a %in% b ==T){
#       levels[i] <- j 
#       break
#     }
#   }
#   
# }
# test8_km_mini$levels <- levels
# 
# test8_km_mini$druggable <- 0
# test8_km_mini$druggable[test8_km_mini$gene %in% druggable$V1] <- 1
# test8_km_mini$cysteine <- 0
# test8_km_mini$cysteine[test8_km_mini$gene %in% cysteine$GENES] <- 1
# 
# 
# 
# 
# test8_hit_mini <- toppNetHITS(mini_logfc,mini_network)
# test8_hit_reac <- toppNetHITS(persudo_logfc,reac_network)
# 
# 
# #9. hidden
# test9_mini <- hiddenNodes(mini_logfc,mini_network,alpha=1,direction = "upstream")
# 


#######load matt's stat model
setwd("~/GITCLONES/stop-sle-statmodel-clustering/combined")
#run install-from-this-source.sh first to get bmtongs package
source("init.R")
all.responses.table = U.R()
responses.table = all.responses.table[comparison=="LNvHC"]
q_filter_gene <- responses.table[responses.table$q <0.2,]

setwd("~/GITCLONES/stop-sle-expression-matrix/")
source(sprintf("%s/stop-sle-expression-matrix/load-rsem.R", Sys.getenv("GITCLONES")))
set.seed(1)
convert <- names(stopsle.rtg)
names(convert) <- stopsle.rtgtest8_km_mini$gene <- unname(mini_to_reac1[test8_km_mini$Node_identifier])

levels <- numeric()
n <- dim(test8_km_mini)[1]
levels <- rep(0,as.numeric(n))
for (i in 1:length(levels)){
  a <-test8_km_mini$gene[i]
  for(j in 1:10){
    b <- up_list[[j]]
    if( a %in% b ==T){
      levels[i] <- j 
      break
    }
  }
  
}
test8_km_mini$levels <- levels

test8_km_mini$druggable <- 0
test8_km_mini$druggable[test8_km_mini$gene %in% druggable$V1] <- 1
test8_km_mini$cysteine <- 0
test8_km_mini$cysteine[test8_km_mini$gene %in% cysteine$GENES] <- 1





new_symbol <- unname(stopsle.rtg[as.character(q_filter_gene$ensembl)])
q_filter_gene$new.symbol <- new_symbol
QGeneList <- split( q_filter_gene, f = q_filter_gene$cell.type)

#####################updatream
level_list <- downstream(c("KLF2","KLF4"),mini_network,how.far.up = 10)
up_list <- list()
net <- readnet(level_list)
for(i in 1:length(net)){

 a <- list(as.character(unname(net[[i]][,1])))
 names(a) <- paste("level",i,sep = "_")
  up_list <- append(up_list,a)
  
}
#3.
test3_mini$gene <- unname(mini_to_reac1[test3_mini$Node_identifier])

levels <- numeric()
n <- dim(test3_mini)[1]
levels <- rep(0,as.numeric(n))
for (i in 1:length(levels)){
  a <-test3_mini$gene[i]
  for(j in 1:10){
    b <- up_list[[j]]
  if( a %in% b ==T){
   levels[i] <- j 
   break
  }
  }

}
test3_mini$levels <- levels

#library(data.table)
#druggable <- fread("druggable.txt",header = F)
#cysteine <- fread("cysteine_cleavage_potential_target_list.csv")
test8_km_mini$gene <- unname(mini_to_reac1[test8_km_mini$Node_identifier])

levels <- numeric()
n <- dim(test8_km_mini)[1]
levels <- rep(0,as.numeric(n))
for (i in 1:length(levels)){
  a <-test8_km_mini$gene[i]
  for(j in 1:10){
    b <- up_list[[j]]
    if( a %in% b ==T){
      levels[i] <- j 
      break
    }
  }
  
}
test8_km_mini$levels <- levels

test8_km_mini$druggable <- 0
test8_km_mini$druggable[test8_km_mini$gene %in% druggable$V1] <- 1
test8_km_mini$cysteine <- 0
test8_km_mini$cysteine[test8_km_mini$gene %in% cysteine$GENES] <- 1


test8_km_mini$gene <- unname(mini_to_reac1[test8_km_mini$Node_identifier])

levels <- numeric()
n <- dim(test8_km_mini)[1]
levels <- rep(0,as.numeric(n))
for (i in 1:length(levels)){
  a <-test8_km_mini$gene[i]
  for(j in 1:10){
    b <- up_list[[j]]
    if( a %in% b ==T){
      levels[i] <- j 
      break
    }
  }
  
}
test8_km_mini$levels <- levels

test8_km_mini$druggable <- 0
test8_km_mini$druggable[test8_km_mini$gene %in% druggable$V1] <- 1
test8_km_mini$cysteine <- 0
test8_km_mini$cysteine[test8_km_mini$gene %in% cysteine$GENES] <- 1



test3_mini$druggable <- 0
test3_mini$druggable[test3_mini$gene %in% druggable$V1] <- 1
test3_mini$cysteine <- 0
test3_mini$cysteine[test3_mini$gene %in% cysteine$GENES] <- 1




#4.
test4_mini$gene <- unname(mini_to_reac1[test4_mini$Node_identifier])

levels <- numeric()
n <- dim(test4_mini)[1]
levels <- rep(0,as.numeric(n))
for (i in 1:length(levels)){
  a <-test4_mini$gene[i]
  for(j in 1:10){
    b <- up_list[[j]]
    if( a %in% b ==T){
      levels[i] <- j 
      break
    }
  }table(test9_mini$druggable[1:100])
  
}
test4_mini$levels <- levels

test4_mini$druggable <- 0
test4_mini$druggable[test4_mini$gene %in% druggable$V1] <- 1
test4_mini$cysteine <- 0
test4_mini$cysteine[test4_mini$gene %in% cysteine$GENES] <- 1


#5.
test5_mini$gene <- unname(mini_to_reac1[test5_mini$Node_identifier])

levels <- numeric()
n <- dim(test5_mini)[1]
levels <- rep(0,as.numeric(n))
for (i in 1:length(levels)){
  a <-test5_mini$gene[i]
  for(j in 1:10){
    b <- up_list[[j]]
    if( a %in% b ==T){
      levels[i] <- j 
      break
    }
  }
  
}
test5_mini$levels <- levels

test5_mini$druggable <- 0
test5_mini$druggable[test5_mini$gene %in% druggable$V1] <- 1
test5_mini$cysteine <- 0
test5_mini$cysteine[test5_mini$gene %in% cysteine$GENES] <- 1



#6.

test6_mini$gene <- unname(mini_to_reac1[test6_mini$Node_identifier])

levels <- numeric()
n <- dim(test6_mini)[1]
levels <- rep(0,as.numeric(n))
for (i in 1:length(levels)){
  a <-test6_mini$gene[i]
  for(j in 1:10){
    b <- up_list[[j]]
    if( a %in% b ==T){
      levels[i] <- j 
      break
    }
  }
  
}
test6_mini$levels <- levels

test6_mini$druggable <- 0
test6_mini$druggable[test6_mini$gene %in% druggable$V1] <- 1
test6_mini$cysteine <- 0
test6_mini$cysteine[test6_mini$gene %in% cysteine$GENES] <- 1



#rev
test6_mini_rev$gene <- unname(mini_to_reac1[test6_mini_rev$Node_identifier])

levels <- numeric()
n <- dim(test6_mini_rev)[1]
levels <- rep(0,as.numeric(n))
for (i in 1:length(levels)){
  a <-test6_mini_rev$gene[i]
  for(j in 1:10){
    b <- up_list[[j]]
    if( a %in% b ==T){
      levels[i] <- j 
      break
    }
  }
  
}
test6_mini_rev$levels <- levels

test6_mini_rev$druggable <- 0
test6_mini_rev$druggable[test6_mini_rev$gene %in% druggable$V1] <- 1
test6_mini_rev$cysteine <- 0
test6_mini_rev$cysteine[test6_mini_rev$gene %in% cysteine$GENES] <- 1

table(test6_mini_rev$druggable[1:100])




#7.

test7_mini$gene <- unname(mini_to_reac1[test7_mini$Node_identifier])

levels <- numeric()
n <- dim(test7_mini)[1]
levels <- rep(0,as.numeric(n))
for (i in 1:length(levels)){
  a <-test7_mini$gene[i]
  for(j in 1:10){
    b <- up_list[[j]]
    if( a %in% b ==T){
      levels[i] <- j 
      break
    }
  }
  
}
test7_mini$levels <- levels

test7_mini$druggable <- 0
test7_mini$druggable[test7_mini$gene %in% druggable$V1] <- 1
test7_mini$cysteine <- 0
test7_mini$cysteine[test7_mini$gene %in% cysteine$GENES] <- 1


#rev
test7_mini_rev$gene <- unname(mini_to_reac1[test7_mini_rev$Node_identifier])

levels <- numeric()
n <- dim(test7_mini_rev)[1]
levels <- rep(0,as.numeric(n))
for (i in 1:length(levels)){
  a <-test7_mini_rev$gene[i]
  for(j in 1:10){
    b <- up_list[[j]]
    if( a %in% b ==T){
      levels[i] <- j 
      break
    }
  }
  
}
test7_mini_rev$levels <- levels

test7_mini_rev$druggable <- 0
test7_mini_rev$druggable[test7_mini_rev$gene %in% druggable$V1] <- 1
test7_mini_rev$cysteine <- 0
test7_mini_rev$cysteine[test7_mini_rev$gene %in% cysteine$GENES] <- 1
table(test7_mini_rev$druggable[1:100])

#8.
test8_km_mini$gene <- unname(mini_to_reac1[test8_km_mini$Node_identifier])

levels <- numeric()
n <- dim(test8_km_mini)[1]
levels <- rep(0,as.numeric(n))
for (i in 1:length(levels)){
  a <-test8_km_mini$gene[i]
  for(j in 1:10){
    b <- up_list[[j]]
    if( a %in% b ==T){
      levels[i] <- j 
      break
    }
  }
  
}
test8_km_mini$levels <- levels

test8_km_mini$druggable <- 0
test8_km_mini$druggable[test8_km_mini$gene %in% druggable$V1] <- 1
test8_km_mini$cysteine <- 0
test8_km_mini$cysteine[test8_km_mini$gene %in% cysteine$GENES] <- 1


test8_hit_mini$gene <- unname(mini_to_reac1[test8_hit_mini$Node_identifier])

levels <- numeric()
n <- dim(test8_hit_mini)[1]
levels <- rep(0,as.numeric(n))
for (i in 1:length(levels)){
  a <-test8_hit_mini$gene[i]
  for(j in 1:10){
    b <- up_list[[j]]
    if( a %in% b ==T){
      levels[i] <- j 
      break
    }
  }
  
}
test8_hit_mini$levels <- levels

test8_hit_mini$druggable <- 0
test8_hit_mini$druggable[test8_hit_mini$gene %in% druggable$V1] <- 1
test8_hit_mini$cysteine <- 0
test8_hit_mini$cysteine[test8_hit_mini$gene %in% cysteine$GENES] <- 1




#9.

test9_mini$gene <- unname(mini_to_reac1[test9_mini$node])

levels <- numeric()
n <- dim(test9_mini)[1]
levels <- rep(0,as.numeric(n))
for (i in 1:length(levels)){
  a <-test9_mini$gene[i]
  for(j in 1:10){
    b <- up_list[[j]]
    if( a %in% b ==T){
      levels[i] <- j 
      break
    }
  }
  
}
test9_mini$levels <- levels

test9_mini$druggable <- 0
test9_mini$druggable[test9_mini$gene %in% druggable$V1] <- 1
test9_mini$cysteine <- 0
test9_mini$cysteine[test9_mini$gene %in% cysteine$GENES] <- 1
table(test9_mini$druggable[1:100])






################################################




########################################################################################################
#Step 4: druggable test, and cysteine test
# Use data table provided here and add a column to identify if the gene is druggable and potential target
######################################################################

library(data.table)
druggable <- fread("druggable.txt",header = F)
cysteine <- na.omit(data.frame(fread("cysteine_cleavage_potential_target_list.csv")))

test_mini$druggable <- 0
test_mini$druggable[test_mini$gene %in% druggable$V1] <- 1
test_mini$cysteine <- 0
test_mini$cysteine[test_mini$gene %in% cysteine$GENES] <- 1
test_mini$source <- "minibase"

# 
# test_reac$druggable <- 0
# test_reac$druggable[test_reac$hypothesis %in% druggable$V1] <- 1
# test_reac$cysteine <- 0
# test_reac$cysteine[test_reac$hypothesis %in% cysteine$GENES] <- 1
# test_reac$source <- "reactome"

test_nwpg_mini <- test6_mini


mini_sign <- test_mini[test_mini$final_rank <= 500,]
mini_sign <- mini_sign[mini_sign$druggable ==1 | mini_sign$cysteine == 1,]
mini_sign <- na.omit(mini_sign$hypothesis)

mini_nwpg <- test_nwpg_mini[1:500,]
mini_nwpg <- mini_nwpg[mini_nwpg$druggable ==1 | mini_nwpg$cysteine == 1,]
mini_nwpg <- na.omit(mini_nwpg$Node_identifier)

sources <- c("-2007163176","-2066067575")

sign_TieDie <- TieDie(targets = sources,sources= mini_sign,mini_network,directed = T)
nwpg_TieDie <- TieDie(targets =sources, sources =  mini_nwpg,mini_network,directed=T)

sub_sign <- sign_TieDie[[1]]
sub_sign$network_object_id1 <-  unname(mini_to_reac1[sub_sign$network_object_id1])
sub_sign$network_object_id2 <- unname(mini_to_reac1[sub_sign$network_object_id2])

sub_nwpg <- nwpg_TieDie[[1]]
sub_nwpg$network_object_id1 <- unname(mini_to_reac1[sub_nwpg$network_object_id1])
sub_nwpg$network_object_id2 <- unname(mini_to_reac1[sub_nwpg$network_object_id2])


gene_mininetwork <- mini_network
gene_mininetwork <- readnet(mini_network)

test <- data.frame(gene= c("KLF2","KLF4"),
                   logfc = c(10,10))
view.subnetwork(sub_nwpg,startnodes = test,directed = T,edgeArrow="triangle")
view.subnetwork(sub_sign,startnodes = test,directed = T,edgeArrow="triangle")


save(test_mini,test_reac,test_nwpg_mini,file = "klf2_klf4_upstream_testing.RData")

View(test_mini[test_mini$druggable ==1 & test_mini$cysteine ==1,])












############################################filter by translational factors
mini_network_tf <- network[-grep("transcription",network$mechanism_name,ignore.case = T),]
test_mini <- SigNet(mini_logfc,mini_network_tf)
test_mini$gene <- unname(mini_to_reac1[test_mini$hypothesis])

test_nwpg_mini <- networkPropagation(mini_logfc,mini_network_tf,directed = T)
test_nwpg_mini$gene <- unname(mini_to_reac1[test_nwpg_mini$Node_identifier])


test_mini$druggable <- 0
test_mini$druggable[test_mini$gene %in% druggable$V1] <- 1
test_mini$cysteine <- 0
test_mini$cysteine[test_mini$gene %in% cysteine$GENES] <- 1

test_nwpg_mini$druggable <- 0
test_nwpg_mini$druggable[test_nwpg_mini$gene %in% druggable$V1] <- 1
test_nwpg_mini$cysteine <- 0
test_nwpg_mini$cysteine[test_nwpg_mini$gene %in% cysteine$GENES] <- 1

mini_sign <- test_mini[test_mini$final_rank <= 1000,]
mini_sign <- mini_sign[mini_sign$druggable ==1 | mini_sign$cysteine == 1,]
mini_sign <- na.omit(mini_sign$hypothesis)

mini_nwpg <- test_nwpg_mini[1:1000,]
mini_nwpg <- mini_nwpg[mini_nwpg$druggable ==1 | mini_nwpg$cysteine == 1,]
mini_nwpg <- na.omit(mini_nwpg$Node_identifier)

sources <- c("-2007163176","-2066067575")

sign_TieDie <- TieDie(targets = sources,sources= mini_sign,mini_network_tf,directed = T)
nwpg_TieDie <- TieDie(targets =sources, sources =  mini_nwpg,mini_network_tf,directed=T)

sub_sign <- sign_TieDie[[1]]
sub_sign$network_object_id1 <-  unname(mini_to_reac1[sub_sign$network_object_id1])
sub_sign$network_object_id2 <- unname(mini_to_reac1[sub_sign$network_object_id2])

sub_nwpg <- nwpg_TieDie[[1]]
sub_nwpg$network_object_id1 <- unname(mini_to_reac1[sub_nwpg$network_object_id1])
sub_nwpg$network_object_id2 <- unname(mini_to_reac1[sub_nwpg$network_object_id2])


test <- data.frame(gene= c("KLF2","KLF4"),
                   logfc = c(10,10))
view.subnetwork(sub_nwpg,startnodes = test,directed = T,edgeArrow="triangle")
view.subnetwork(sub_sign,startnodes = test,directed = T,edgeArrow="triangle")






#############################################################wrting function for upstream regulating





a <- mini_network
a$gene1 <- unname(mini_to_reac1[a$network_object_id1])
a$gene2 <- unname(mini_to_reac1[a$network_object_id2])
a <- na.omit(a[a$edgetype !=0,])


lvl1 <- a[a$gene2 == "KLF2",]
a <- lvl1
for (i in 1:4){
 
  a <- klf2[klf2$gene2 %in% a$gene1,]
  assign(paste("lvl",i+1,sep = ""),a)
  
}


k4lvl1 <- a[a$gene2 == "KLF4",]


a <- k4_lvl1
for (i in 1:4){
  
  a <- klf2[klf2$gene2 %in% a$gene1,]
  assign(paste("k4lvl",i+1,sep = ""),a)
  
}

k2gene <- as.character(c(lvl1$gene1,lvl2$gene1,lvl3$gene1,lvl4$gene1,lvl5$gene1))
k4gene <- as.character(c(k4lvl1$gene1,k4lvl2$gene1,k4lvl3$gene1,k4lvl4$gene1,k4lvl5$gene1))
k2gene <- unique(k2gene[k2gene %in% k4gene])
k4gene <- unique(k4gene[k4gene %in% k2gene])
common <- k2gene
common <- append(common,c("KLF2","KLF4"))

#process again
a <- mini_network
a$gene1 <- unname(mini_to_reac1[a$network_object_id1])
a$gene2 <- unname(mini_to_reac1[a$network_object_id2])
a <- na.omit(a[a$edgetype !=0,])
a <- a[a$gene1 %in% common | a$gene2 %in% common,]



lvl1 <- a[a$gene2 == "KLF2",]
a <- lvl1
for (i in 1:4){
  
  a <- klf2[klf2$gene2 %in% a$gene1,]
  assign(paste("lvl",i+1,sep = ""),a)
  
}


k4lvl1 <- a[a$gene2 == "KLF4",]
a <- k4lvl1
for (i in 1:4){
  
  a <- klf2[klf2$gene2 %in% a$gene1,]
  assign(paste("k4lvl",i+1,sep = ""),a)
  
}





ll1 <- lvl1
names <- ll1$edgetype
ll1 <- as.character(ll1$gene1)
names(ll1) <- names
k2l2 <- as.character()
for (i in 1:nrow(lvl1)){
  i=2
    ll2 <-  lvl2[lvl2$gene2 == lvl1$gene1[i],]
    names <- ll2$edgetype
    ll2 <- as.character(paste(ll2$gene1,lvl1$gene1[i],sep = "_"))
   names(ll2) <- names
   k2l2 <- append(k2l2,ll2)
   #lvl3 
   
  
  
  
}















l1 <- paste(lvl1$gene1,lvl1$gene2,sep = "_")
names(l1) <- lvl1$edgetype
l2 <- paste(lvl2$gene1,lvl2$gene2,sep = "_")
names(l2) <- lvl2$edgetype

l3 <- paste(lvl3$gene1,lvl3$gene2,sep = "_")
names(l3) <- lvl3$edgetype




















k4l1 <- paste(k4lvl1$gene1,k4lvl1$gene2,sep = "_")
names(k4l1) <- k4lvl1$edgetype
k4l2 <- paste(k4lvl2$gene1,k4lvl2$gene2,sep = "_")
names(k4l2) <- k4lvl2$edgetype

k4l3 <- paste(k4lvl3$gene1,k4lvl3$gene2,sep = "_")
names(k4l3) <- k4lvl3$edgetype



# 1 is what we want

ii=list()
gg=list()
for(i in 1:length(l2)){
i=1
gene1 <- gsub('.*\\_',"",l2[i])
gene2 <- l1[grep(gene1,gsub('\\_.*',"",l1))]
m <- (as.numeric(names(gene1)))*(as.numeric(names(gene2)))

# k4gene1 <- gsub('.*\\_',"",k4l2[i])
# k4gene2 <- k4l1[grep(k4gene1,gsub('\\_.*',"",k4l1))]
# n <- (as.numeric(names(k4gene1)))*(as.numeric(names(k4gene2)))
#         
if(m == 1){
  print(paste(l2[i],gene2,"is a up regulating pathway for KLF2 and KLF4"))
ii <- append(ii,list(l2[i],gene2))

gene3 <- gsub('\\_.*',"",l2[i])
for (j in 1: length(l3))
gene4 <- l3[grep(gene3,gsub('.*\\_',"",l3))]
n <- as.numeric(names(gene4))*m
if(n ==1){
  print(paste(l3[j],gene4,l2[i],gene2,"is a up regulating pathway for KLF2 and KLF4"))
  gg <- append(gg,list(l3[j],gene4))
  
}
}
  
}
}
#gene2 <- gsub('.*\\_',"",l3[i])
#gene3 <- l2[grep(gene2,gsub('\\_.*',"",l2))]
#n <- (as.numeric(names(gene2)))*(as.numeric(names(gene3)))














b <- sub_nwpg[sub_nwpg$network_object_id2 == "KLF4"&sub_nwpg$edgetype !=0,]
