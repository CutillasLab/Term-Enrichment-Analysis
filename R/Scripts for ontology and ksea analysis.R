

## Functions required for ontology and kinase substrate enrichment analysis


########## Dot Plot of selected pathway s#########################################
plots.selected.ont.db <- function(mypathway){
  temp.phospho <- merge(df.sen.phospho,df.res.phospho,by=c("pathway.db","drug"),all = T)
  temp.phospho$enrichment.x[is.na(temp.phospho$enrichment.x)] <- 1
  temp.phospho$enrichment.y[is.na(temp.phospho$enrichment.y)] <- 1
  temp.phospho$delta.enrich <- log2(temp.phospho$enrichment.x)-log2(temp.phospho$enrichment.y)
  temp.phospho <- merge(temp.phospho,df.drug.info,by="drug")
  
  pp1 <- plot.enrichment(target.pathway = mypathway,
                         df.sen = temp.phospho,
                         db = "pdts",
                         input="Phosphoproteomics AML")
  pp1
  #Proteomics
  
  temp.prot <- merge(df.sen.protein,df.res.protein,by=c("pathway.db","drug"),all = T)
  temp.prot$enrichment.x[is.na(temp.prot$enrichment.x)] <- 1
  temp.prot$enrichment.y[is.na(temp.prot$enrichment.y)] <- 1
  temp.prot$delta.enrich <- log2(temp.prot$enrichment.x)-log2(temp.prot$enrichment.y)
  temp.prot <- merge(temp.prot,df.drug.info,by="drug")
  
  #temp.prot <- subset(temp.prot,temp.prot$FDR.x<0.05 | temp.prot$FDR.y<0.05)
  
  levels(as.factor(temp.prot$ont.database.x))
  
  pp2 <- plot.enrichment(target.pathway = mypathway,
                         df.sen = temp.prot,
                         db = "nci",
                         input="proteomics AML")
  pp2
  #Transcriptomics
  
  
  temp.rna <- merge(df.sen.rna,df.res.rna,by=c("pathway.db","drug"),all = T)
  temp.rna$enrichment.x[is.na(temp.rna$enrichment.x)] <- 1
  temp.rna$enrichment.y[is.na(temp.rna$enrichment.y)] <- 1
  temp.rna$delta.enrich <- log2(temp.rna$enrichment.x)-log2(temp.rna$enrichment.y)
  temp.rna <- merge(temp.rna,df.drug.info,by="drug")
  
  #temp.rna <- subset(temp.rna,temp.rna$FDR.x<0.05 | temp.rna$FDR.y<0.05)
  
  levels(as.factor(temp.rna$ont.database.x))
  
  pp3 <- plot.enrichment(target.pathway = mypathway,
                         df.sen = temp.rna,
                         db = "tf.targets.omnipath",
                         input="Transcriptomics AML")
  
  pp3
  
  return(list(pp1,pp2,pp3))
}
########## Plot enrichment #######################################################
# 
plot.enrichment <- function(target.pathway,
                            df.sen,
                            db="hallmark.genes",
                            input=""){
  df.sen$delta.enrich <- scale(df.sen$delta.enrich)
  xsen <- df.sen[df.sen$Target.pathway==target.pathway &
                   df.sen$ont.database.x==db,]
  
  xsen.cast <- (dcast(xsen, pathway.db~drug.target,
                      value.var =  "delta.enrich" ,
                      fun.aggregate = mean ))
  xsen.cast[is.na(xsen.cast)] <- 0
  ss <- apply(xsen.cast[,2:ncol(xsen.cast)],1,function(x){sum(x!=0)})
  xsen.cast <- subset(xsen.cast,xsen.cast$pathway.db!=0
                      & ss>2)
  rownames(xsen.cast) <- xsen.cast$pathway.db
  df.pp <- as.matrix(xsen.cast[,2:ncol(xsen.cast)])
  df.pp <- df.pp[,which( colnames(df.pp)!="NA")]
  
  palette.breaks<- seq(0,7,0.1)
  color.palette <- colorRampPalette(c('white','orange','red'))(length(palette.breaks) - 1)
  
  hm <- gplots::heatmap.2((df.pp))
  
  plot.sen <- ggplot(xsen,
                     aes(x=pathway.db,y=drug.target))+
    geom_point(aes(size=abs(delta.enrich), color=delta.enrich))+
    #geom_jitter(aes(color=Target.pathway),width = 0.2)+
    #scale_color_manual(values = c23)+
    scale_x_discrete(limits=rownames(df.pp)[hm$rowInd])+
    scale_y_discrete(limits=colnames(df.pp)[hm$colInd])+
    scale_color_gradient2(low="seagreen",mid="white",high = "red",midpoint = 0)+
    theme_bw()+
    theme(legend.position = "left",
          axis.text.x = element_text(angle = 90, hjust=1))+
    labs(title = paste("Enrichment of markers for drugs against",target.pathway,db ),
         subtitle = paste("Input =",input))
  plot.sen
  return(plot.sen)
}
########## Dendogram of delta enrichment values #################################################################
plot.dendo <- function(df.s,input.data="",plot.legend=TRUE){
  df.sim.scores <- data.frame(dcast(df.s,drug.target ~ pathway.db ,value.var = "delta.enrich", mean))
  df.sim.scores <- subset(df.sim.scores,is.na(df.sim.scores$drug.target)==F)
  rownames(df.sim.scores) <- df.sim.scores$drug.target
  df.sim.scores[is.na(df.sim.scores)] <- 0
  df.sim.scores <-  df.sim.scores[,2:ncol(df.sim.scores)]

  drugs <- rownames(df.sim.scores)
  mycols <- character()
  shapes <- character()
  r <- 1
  for (drug in drugs) {
    mycols[r]<- df.drug.info[df.drug.info$drug.target==drug,"mycols"]# $Target.pathway[r]
    shapes[r] <- df.drug.info[df.drug.info$drug.target==drug,"shapes"]
    r <- r+1
  }
  dd <- dist(((df.sim.scores)), method = "euclidean")
  hc <- hclust(dd, method = "median")
  
  plot(as.phylo(hc),   type = "fan",     tip.color = mycols,
       cex = .7,no.margin = TRUE,       main=input.data)
  if (plot.legend==TRUE){
    legend( "topright",legend = unique(df.cols$pathways), col = unique((df.cols$cols)),lty= 1,lwd = 5,cex=.7)
  }
  
}

########## Drug Similarity Score calculation #############################################################################################
# Calculate similarity scores as the pearson r value of delta enrichment values between drugs
drug.similarity.scores <- function(df.s,input.data=""){
  df.sim.scores <- data.frame(dcast(df.s,drug.target ~ pathway.db ,value.var = "delta.enrich", mean))
  df.sim.scores <- subset(df.sim.scores,is.na(df.sim.scores$drug.target)==F)
  rownames(df.sim.scores) <- df.sim.scores$drug.target
  df.sim.scores[is.na(df.sim.scores)] <- 0
  df.sim.scores <-  df.sim.scores[,2:ncol(df.sim.scores)]
  df.cor <- data.frame(cor(t(df.sim.scores[,2:ncol(df.sim.scores)])))
  df.cor <- melt(as.matrix(df.cor))
  colnames(df.cor) <- c("drug1","drug.2","similarity")
  return(df.cor)
}
########## Plot similarity Scores ##################################
# Plot similarity scores
plot.similarity.score.for.drug <- function(drug){
  df.s <- df.sim[grep(drug,df.sim$drug1),]
  df.s1 <- df.s[order(-df.s$similarity),]
  df.s2 <- df.s[order(df.s$similarity),]
  #df.ss <- rbind.data.frame(  df.s1[1:10,],df.s2[1:10,])
  #head(df.ss[,1:10],n=20)
  
  ps3 <- ggplot(df.s1[1:15,],aes(x=drug.2,y=similarity))+
    geom_bar(stat = "identity")+
    #scale_fill_manual(values=c23)+
    #scale_x_discrete(limits=c(df.s2$drug.target[1:10],rev(df.s1$drug.target[1:10])))+
    scale_x_discrete(limits=rev(df.s1$drug.2[1:15]))+
    ylim(c(0,1))+
    labs(title = paste("Similarity score for",drug))+
    theme_classic()+coord_flip()
  return(ps3)
}

plot.similarity.score.for.drug.2 <- function(drug){
  df.s <- df.sim[grep(drug,df.sim$drug1),]
  df.s$scaled.similarity <- scale(df.s$similarity)
  
  
  mm <- median(df.s$scaled.similarity)
  ss <- sd(df.s$scaled.similarity)
  
  #df.ss <- rbind.data.frame(  df.s1[1:10,],df.s2[1:10,])
  #head(df.ss[,1:10],n=20)
  df.s$pvalue <- 1
  for (i in 1:nrow(df.s)){
  res <- t.test(df.s$scaled.similarity, mu = df.s$scaled.similarity[i])
  df.s$pvalue[i] <- res$p.value
  }
  df.s1 <- df.s[order(-df.s$similarity),]
  df.s2 <- df.s[order(df.s$similarity),]
  ps3 <- ggplot(df.s,aes(x=scaled.similarity,y=drug.2))+
    geom_point(aes(size=-log10(pvalue), color=-log10(pvalue)))+
    geom_text_repel(data = df.s1[1:20,], 
                    aes(x=scaled.similarity,y=drug.2, label=drug.2,
                        color=-log10(pvalue)),
                    size=2.5)+
    scale_color_gradient(high = "purple4",low="grey")+
    
    labs(title = paste("Similarity score for",drug))+
    theme_classic()+
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())+
    geom_vline(xintercept = 0, linetype=2)+
    geom_vline(xintercept = 1, linetype=3)+
    geom_vline(xintercept = -1, linetype=3)
  
  ps3
  return(ps3)
}

########## Onoloty and KS enrichment analysis from list of EMDRs ###############################################################

Enrichment.from.list <- function(list.of.peptides,background.list,prot_db=c("edges","pdts","psite","signor","reactome","process","function","location", "kegg","nci","tf.targets.omnipath","hallmark.genes"),
                                 is.ksea=FALSE){

  # Returns enrichment of phosphosites, proteins or transcripts in pathways ontologies or kinase-phosphosite relationsips
  #
  # list of peptides is the list of proteins, phosphosites or transcripts to be searched (transcripts needt o be converted to Uniprot names)
  # Accepts Uniprot names or phosphosite names in the form gene(xy), were gene = uniprot gene name; x=S, T or Y; y= amino acid position of phosphorylation
  # background.list is the background list of of proteins, phosphosites or transcripts
  # prot_db = any suitable ontology or kinase-substrate database
  # is.ksea = change to TRUE for analysing phosphoproteomics data against kinase-substrates
  # draw.tables.and.plots = if TRUE draws volcano plots of the enrichment data, requires ggplot2
  #
  # Enrichment = (q/k)/(j/m), where:
  #   q = peptides/proteins in list.of.peptides with a match in the ontology/K-S relationship dataset
  #   k = number of peptides/proteins in list.of.peptides
  #   j = peptides/proteins in background.list with a match in the ontology/K-S relationship dataset
  #   m = number of peptides/proteins in background.list
  #
  # p-values are calculated using the hypergeometric function and adjusted by FDR method
  #



  if (is.ksea==FALSE){
    x1 <-  grepl(")", as.character(list.of.peptides),fixed=T)
    x2 <-  grepl("(", as.character(list.of.peptides),fixed = T)

    if (x1[1]==TRUE & x2[1]==TRUE){
        phosphosite.names <- list.of.peptides
        list.of.peptides <- accessions.from.phosphosite.names(phosphosite.names)
      }
      x3 <-  grepl("..", as.character(list.of.peptides),fixed = T)
    if (x3[1]==TRUE){
        phosphosite.names <- list.of.peptides
        list.of.peptides <- accessions.from.phosphosite.names(phosphosite.names)
      }
      x1 <-  grepl(")", as.character(background.list),fixed=T)
      x2 <-  grepl("(", as.character(background.list),fixed = T)
    if (x1[1]==TRUE & x2[1]==TRUE){
        background.list <- accessions.from.phosphosite.names(background.list)
      }
      x3 <-  grepl("..", as.character(background.list),fixed = T)
    if (x3[1]==TRUE){
        background.list <- accessions.from.phosphosite.names(background.list)
      }
  }

  df.ks <- get.protein.set.data(prot_db)
  df.ks <- df.ks[order(-df.ks[,2]),]
  nr <- nrow(df.ks)

  pathway <- character(nr)
  pvalues <- numeric(nr)
  enrichment <- numeric(nr)
  counts <- numeric(nr)
  counts.bg <- numeric(nr)
  FDR <- numeric(nr)
  proteins <- character(nr)
  
  m <- length(background.list)
  n <- length(list.of.peptides)
  df.peptides <- data.frame(peptides=list.of.peptides,
                            proteins=list.of.peptides)
  bg.size <- numeric(nr)
  data.size <- numeric(nr)
  r=1
  for (r in 1:nr) {
    mym <- df.ks[r,2]
    kinase <- df.ks[r,1]
    if (is.na(mym) == F){
      if(mym>2){
        substrates <- as.character(df.ks[r,3])
        ss <- c(unlist(strsplit(substrates,";")))
        if (is.ksea==TRUE | is.ksea==T){
          ss <- paste(ss,";",sep = "")
        }
        start.time <- Sys.time()
        k <- n#length(ss)
        q <- length(intersect(ss,list.of.peptides))
        j <- length(intersect(ss,background.list))
        if (q>0 & j>0){
          prots <- intersect(ss,list.of.peptides)
          peptides <- df.peptides[df.peptides$proteins %in% prots,"peptides"]
          pvalue <- 1-phyper(q-1,j,m-j,k,lower.tail = TRUE, log.p = F)
          pathway[r] <- as.character(kinase)
          pvalues[r] <- pvalue
          enrichment[r] <- round((q/k)/(j/m), digits = 2)
          counts[r] <- q
          data.size[r] <- k
          counts.bg[r] <- j
          bg.size[r] <- m
          proteins[r] <- paste(peptides,collapse = ";")
          
          }
      }
      r=r+1
    }
  }
  results <- data.frame(pathway, pvalues,enrichment,counts, data.size,counts.bg,bg.size, proteins)
  results <- na.omit(results)
  results <- results[order(-results$enrichment),]
  results <- subset(results, pvalues!=0 & counts>1)
  results$FDR <- p.adjust(results$pvalues, method = "fdr")

  return(results)
}

get.protein.set.data <- function(dataset){
  ## data sets
  edges <- "https://www.dropbox.com/s/ttmzd40mnjgh1iu/edges.csv?dl=1"
  pdts <- "https://www.dropbox.com/s/86jfnayv0qa1n2q/pdts.csv?dl=1"
  psite <- "https://www.dropbox.com/s/eb1qoofz793f4tq/psite.csv?dl=1"
  signor <- "https://www.dropbox.com/s/alpbq880emz1z2t/signor.csv?dl=1"
  reactome <- "https://www.dropbox.com/s/jdcc1355cz73mmi/reactome.csv?dl=1"
  process <- "https://www.dropbox.com/s/z8zef96q45vi0je/process.csv?dl=1"
  myfunction <- "https://www.dropbox.com/s/ev50zpd0g41pds3/function.csv?dl=1"
  location <- "https://www.dropbox.com/s/js8vwyamlhbnqqf/location.csv?dl=1"
  nci <- "https://www.dropbox.com/s/fe8t4nyhbljsn5y/nci.csv?dl=1"
  tf.targets.omnipath <- "https://www.dropbox.com/s/hg3slk150l7zd0x/TF%20all%20targets%20omnipath.csv?dl=1"
  kegg <- "https://www.dropbox.com/s/gm2821cmxarv7sx/kegg%20pathways.csv?dl=1"
  hallmark.genes <-  "https://www.dropbox.com/s/wqvnoalg2v6ufm8/hallmark%20genes%20v71.csv?dl=1"
  dataset.names <-  c("edges","pdts","psite","signor","reactome","process","function","location","nci", "tf.targets.omnipath","kegg", "hallmark.genes")
  datasets <- c(edges, pdts,psite,signor,reactome,process,myfunction,location,nci,tf.targets.omnipath,kegg, hallmark.genes)
  df.datasets <- data.frame(dataset.names,datasets)
  myfile <- as.character(df.datasets[df.datasets$dataset.names==dataset,2])
  if (length(myfile)==0){
    print (paste(dataset, "not found. Check spelling"))
  }
  df.out <- read.csv(myfile)
  return(df.out)
}
accessions.from.phosphosite.names <- function(phosphosite.names){
  prot.data <- read.csv("https://www.dropbox.com/s/q0l7lqsy7nt7tlf/uniprot_names.csv?dl=1")
  nr <- length(phosphosite.names)
  genes <- lapply(phosphosite.names,
                  function(x) strsplit(as.character(x),"(",fixed = TRUE)[[1]][1])
  df.s1 <- data.frame(gene=unlist(genes),acc="acc",stringsAsFactors = F)
  for (i in 1:nr ){
    acc <- "NA"
    gene <- as.character(df.s1$gene[i])
    if (length(gene)>0){
      acc.1 <- subset(prot.data, prot.data$gene.name==gene)
      acc <- as.character(acc.1$acc[1])
    }
    df.s1$acc[i] <- acc
  }
  print(nr)
  print ("Accession names converted.")
  return(df.s1$acc)

}



########## Colors used in plots ####################
c25 <- c("dodgerblue2", "#E31A1C", # red
         "green4",
         "#6A3D9A", # purple
         "#FF7F00", # orange
         "black",
         "skyblue2",
         "palegreen2",
         "#CAB2D6", # lt purple
         "orangered", # lt orange
         "gray70", "khaki2",
         "maroon", "orchid1",  "blue1", "steelblue4",
         "darkturquoise", "green1", "yellow4", "yellow3",
         "darkorange4", "brown", "gold1",
         "#FB9A99", "deeppink1", "royalblue" # lt pink
)