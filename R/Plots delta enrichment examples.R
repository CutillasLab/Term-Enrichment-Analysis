

# Required libraries
library(grid)
library(gridExtra)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(cowplot)
library(readxl)
library("ape")
library(ggdendro)
library(data.table)

# Required functions are in R file named "Scripts for ontology and ksea analysis"

###########################################################################
###########################################################################
# info drugs
df.drug.info <- read.csv("https://www.dropbox.com/s/7ht3yqdj05hw566/Drug%20info%20from%20Shelleckchem_b.csv?dl=1",row.names=1)
rownames(df.drug.info) <- df.drug.info$drug
rownames(df.drug.info) <- gsub("-",".",rownames(df.drug.info),fixed = T)
rownames(df.drug.info) <- gsub("(",".",rownames(df.drug.info),fixed = T)
rownames(df.drug.info) <- gsub(")",".",rownames(df.drug.info),fixed = T)
###########################################################################

## Load resutls of ontology analysis ### Input is Suppl Dataset 5

f <- "https://www.dropbox.com/s/odbxlv7s5tw8627/Suppl%20Dataset%205%20Systematic%20ontology%20analysis%20of%20EMDRs.xlsx?dl=1"
df.res.phospho <- data.frame(readxl::read_excel(f,"Resistance Phospho AML" ))
df.sen.phospho <- data.frame(readxl::read_excel(f,"Sensitivity Phospho AML" ))
df.res.protein <- data.frame(readxl::read_excel(f,"Resistance Protein AML" ))
df.sen.protein <- data.frame(readxl::read_excel(f,"Sensitivity Protein AML" ))
df.res.rna <- data.frame(readxl::read_excel(f,"Resistance Transcript AML" ))
df.sen.rna <- data.frame(readxl::read_excel(f,"Sensitivity Transcript AML" ))

# Plot example
mypathway <- "PI3K/MTOR signaling"
ppx.pi3k <- plots.selected.ont.db("PI3K/MTOR signaling")
plot.all <- plot_grid(plotlist = ppx.pi3k,nrow = 3,rel_heights = c(1,1.6,1))
ggsave(filename = "Kinases and pathways enrichment PI3K inhibitors.pdf",plot = plot.all,
       width = 10, height = 15)

############################

# Add color to pathway
pathways <- levels(as.factor(df.drug.info$Target.pathway))
df.cols <- data.frame(pathways,cols=c25[1:length(pathways)])
df.s$cols <- "white"
for (r in 1:nrow(df.s)) {
  pp <- df.s$Target.pathway[r]
  df.s$cols[r] <- df.cols[df.cols$pathways==pp,"cols"]
}
df.drug.info$mycols <- c(rep("white", nrow(df.drug.info)))
df.drug.info$shapes <- c(rep(1, nrow(df.drug.info)))
for (r in 1:nrow(df.drug.info)) {
  #drug < df.drug.info$drug.target[r]
  pp <-  df.drug.info$Target.pathway[r]
  df.drug.info$mycols[r] <- df.cols[df.cols==pp,"cols"]
  df.drug.info$shapes[r] <- rownames(df.cols[df.cols==pp,])
}
############################################


# clustering
head(temp.phospho)
temp.phospho <- merge(df.sen.phospho,df.res.phospho,by=c("pathway.db","drug"),all = T)
temp.phospho$enrichment.x[is.na(temp.phospho$enrichment.x)] <- 1
temp.phospho$enrichment.y[is.na(temp.phospho$enrichment.y)] <- 1
temp.phospho$delta.enrich <- log2(temp.phospho$enrichment.x)-log2(temp.phospho$enrichment.y)
temp.phospho <- merge(temp.phospho,df.drug.info,by="drug")
db <- "pdts"

  df.s <- temp.phospho[temp.phospho$Target.pathway !="Other"
                   &temp.phospho$Target.pathway !="Other, kinases"
                   &temp.phospho$ont.database.x==db,]
  #pdf(file = paste("phospho",db,"no legend.pdf"),width = 8, height = 8)
  plot.dendo(df.s=df.s,input.data = paste("Phospho",db), plot.legend = FALSE)
 # dev.off()

################################################################

## Similarity Scores #######################
# Get similarity indeces for all datasets
all.data.in <- c("phospho","protein","rna")
ont.datasets <- list(temp.phospho,temp.prot,temp.rna)
i <- 1
for (data.in in all.data.in ){
  #dat.in <- all.data.in[i]
  ont.dat <- ont.datasets[[i]]
  df.sim.all <- ""
  dbs <- levels(as.factor(ont.dat$ont.database.x))
  for (db in dbs){
    df.s <- ont.dat[ont.dat$ont.database.x==db,]
    df.sim <- drug.similarity.scores(df.s = df.s)
    df.sim$input.data.set <- paste(data.in,db)
    df.sim.all <- rbind.data.frame(df.sim.all,df.sim)
  }
  i <- i+1
  write.csv(df.sim.all,file = paste("Suppl dataset 6 Similarity indeces",dat.in,".csv"))
}

# get similarity indeces for just pdts
db <- "pdts"
df.s <- temp.phospho[temp.phospho$Target.pathway !="Other"
                     &temp.phospho$ont.database.x==db,]
df.sim <- drug.similarity.scores(df.s = df.s) %>% as.data.table

# Plot some examples
pp2 <- plot.similarity.score.for.drug("BYL-719")
pp4 <- plot.similarity.score.for.drug("Rapamycin")
pp8 <- plot.similarity.score.for.drug("GSK1120212")
pp.all <- cowplot::plot_grid(pp2,pp4,pp8,nrow = 1,
                             rel_widths = c(1.3,1.3,1.1))
ggsave(filename = "Similarity socores for selected drugs.pdf",
       plot = pp.all, width = 12, height = 3)

