# Pedro R Cutillas, November 2020

# Example of use of ontology enrichment functions to calcuate ontology enrichment of empirical drug response markers
# shown in Figure 3
# Required functions are in R file named "Scripts for ontology and ksea analysis"

# get background proteins
f.prot <- "https://www.dropbox.com/s/snrjh4psdnbnl2c/ProtQuant.csv?dl=1"
df.prot <- read.csv(f.prot)
bg.proteins <- unique(df.prot$X.1)
# get Empirical Markers of Drug Response
mm <- "https://www.dropbox.com/s/yddjt7e0lrdnyfa/prot_sensitivity_markers_aml.csv?dl=1"
df.markers <- read.csv(mm,row.names=1)
nr <- nrow(df.markers)
ont <- "kegg"
onts <- c("hallmark.genes","process","function","location", "tf.targets.omnipath" ,"nci","reactome","kegg")
df.sen <- ""
df.res <- ""
for (ont in onts){
    for (r in 1:nr){
      # Get markers of drug sensitivity for each drug
      drug <- df.markers$drug[r]
      sen.markers <- df.markers$up.in.sensitive[r]
      m1 <- unlist(strsplit(as.character(sen.markers),"-",fixed = T))
      tryCatch({
        # Enrichment of  sensitivy markers in the specfied ontology 
        ee.table <- Enrichment.from.list(list.of.peptides = m1, background.list = bg.proteins,prot_db = ont)
        en.table$drug <- drug
        en.table$ont.database <- ont
        # merge results for all drugs into a data frame
        df.sen <- rbind.data.frame(en.table,df.sen)
        },error=function(e){} )
    }
    # The same as above but for resistance markers
    for (r in 1:nr){
      # Get markers of drug resistance for each drug
      drug <- df.markers$drug[r]
      res.markers <- df.markers$up.in.resistant[r]
      m1.r <- unlist(strsplit(as.character(res.markers),"-",fixed = T))
      tryCatch({
        ee.table <- Enrichment.from.list(list.of.peptides = m1.r,background.list = bg.proteins,prot_db = ont)
        en.table$drug <- drug
        en.table$ont.database <- ont
        df.res <- rbind.data.frame(en.table,df.res)
      },error=function(e){})
    }
}
write.csv(df.sen,file="Ont enrichment up sensitive.csv")
write.csv(df.res,file="Ont enrichment up resistant.csv")
##################################################################################