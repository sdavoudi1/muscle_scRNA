# Setup ---------------------------

# This is Brendan's code to prepare the DE analysis of the interactome comparing young and aged 
# muscle
library(scales)
library(RColorBrewer)
library(Matrix)
library(scClustViz)
sysDr <- switch(Sys.info()["sysname"],"~/",Windows="D:/")
netDir <- paste0(sysDr,"Dropbox/GDB/CCCnetResources/1705_RuthData/")
inputDataPath <- paste0(sysDr,"Dropbox/GDB/GilbertMuscle/output_v3m_downsampled_perClust2/")
dataPath <- paste0(sysDr,"Dropbox/GDB/GilbertMuscle/output_v4m/")
dataTitle <- "MuscleAging"
dir.create(dataPath,showWarnings=F)


# Load DE results -----------------------
deL <- list()
for (i in grep(paste0(dataTitle,"_selDE_"),list.files(inputDataPath),value=T)) {
  if (grepl("AllaAgedbYoung",i)) { next }
  l <- gsub("MuscleAging_selDE_|.RData","",i)
  load(paste0(inputDataPath,i))
  if (all(sign(new_deMarker[[1]]$`Set aged`$dDR) > 0)) {
    new_deMarker[[1]]$`Set aged`$dDR <- new_deMarker[[1]]$`Set aged`$dDR * -1
    new_deMarker[[1]]$`Set aged`$logGER <- new_deMarker[[1]]$`Set aged`$logGER * -1
  }
  if (all(sign(new_deMarker[[1]]$`Set young`$dDR) < 0)) {
    new_deMarker[[1]]$`Set young`$dDR <- new_deMarker[[1]]$`Set young`$dDR * -1
    new_deMarker[[1]]$`Set young`$logGER <- new_deMarker[[1]]$`Set young`$logGER * -1
  }
  deL[[l]] <- rbind(new_deMarker[[1]]$`Set aged`,new_deMarker[[1]]$`Set young`)
  deL[[l]]$scoreYXoverOX <- -log10(deL[[l]]$qVal) * sign(deL[[l]]$dDR)
  deL[[l]]$scoreYXoverOX[is.na(deL[[l]]$scoreYXoverOX)] <- 0
  deL[[l]]$scoreYXoverOXscaled <- sapply(deL[[l]]$scoreYXoverOX,function(X) 
    X / switch(as.character(X >= 0),
               "TRUE"=max(deL[[l]]$scoreYXoverOX),
               "FALSE"=min(deL[[l]]$scoreYXoverOX) * -1))
  deL[[l]]$scoreYXoverOXscaled[is.na(deL[[l]]$scoreYXoverOXscaled)] <- 0
  if (length(unique(deL[[l]]$scoreYXoverOXscaled)) == 1) { 
    deL[[l]]$scaledFactor <- 6
  } else { 
    deL[[l]]$scaledFactor <- 6
    if (any(deL[[l]]$scoreYXoverOXscaled < 0)) {
      deL[[l]]$scaledFactor[deL[[l]]$scoreYXoverOXscaled < 0] <- 
        cut(deL[[l]]$scoreYXoverOXscaled[deL[[l]]$scoreYXoverOXscaled < 0],5,labels=F)
    }
    if (any(deL[[l]]$scoreYXoverOXscaled > 0)) {
      deL[[l]]$scaledFactor[deL[[l]]$scoreYXoverOXscaled > 0] <- 
        deL[[l]]$scaledFactor[deL[[l]]$scoreYXoverOXscaled > 0] +
        cut(deL[[l]]$scoreYXoverOXscaled[deL[[l]]$scoreYXoverOXscaled > 0],5,labels=F)
    }
  }
}
rm(list=grep("^new",ls(),value=T))
inputGeneList <- unique(Reduce(union,sapply(deL,rownames)))


# prep_ligand_receptor_network -----------------
if (!file.exists(paste0(netDir,"LRNforCCCnetPrediction_editedbyBI_noECM_noNR.RData"))) {
  #### Load interaction data ####
  inxDB_pf1 <- read.table(paste0(netDir,"LigandReceptorInteractions.txt"),header=T,sep="\t",as.is=T)
  for (i in colnames(inxDB_pf1)) {
    inxDB_pf1[which(inxDB_pf1[,i] == "" | inxDB_pf1[,i] == "-"),i] <- NA
  }
  colnames(inxDB_pf1)[colnames(inxDB_pf1) %in% c("AliasA","AliasB")] <- c("nodeA","nodeB")
  
  geneInfo_pf <- read.table(paste0(netDir,"ProteinTypesAttribs_editedbyBI_noECM_noNR.txt"),
                            header=T,sep="\t",as.is=T)
  geneInfo_pf <- geneInfo_pf[!geneInfo_pf$curated_protein_type == "",]
  colnames(geneInfo_pf)[grepl("protein_type",colnames(geneInfo_pf))] <- "protein_type"
  geneInfo_pf <- geneInfo_pf[geneInfo_pf$hgnc_symbol %in% unique(c(inxDB_pf1$nodeA,inxDB_pf1$nodeB)),]
  rownames(geneInfo_pf) <- geneInfo_pf$hgnc_symbol
  
  #### Remove duplicate interactions ####
  inxDB_pf <- as.data.frame(t(apply(inxDB_pf1,1,function(X) {
    temp_order <- unlist(lapply(c(1,3,5,7,12),function(Y) order(X[c("nodeA","nodeB")]) + Y))
    temp_key <- paste(sort(X[c("nodeA","nodeB")]),collapse="_")
    return(c(temp_key,X[c(temp_order[1:8],10:12,temp_order[9:10],15:18)]))
  })),stringsAsFactors=F)
  colnames(inxDB_pf) <- colnames(inxDB_pf1)
  inxDB_pf <- aggregate(inxDB_pf[,2:ncol(inxDB_pf)],by=list(inxDB_pf$key),
                        function(X) paste(unique(unlist(strsplit(X,"|",fixed=T))),collapse="|"))
  colnames(inxDB_pf)[grep("^Group",colnames(inxDB_pf))] <- "key"
  rownames(inxDB_pf) <- inxDB_pf$key
  
  #### Remove self-loops ####
  inxDB_pf <- inxDB_pf[!apply(inxDB_pf,1,function(X) X["nodeA"] == X["nodeB"]),]
  
  #### Remove unannotated proteins ####
  inxDB_pf <- inxDB_pf[inxDB_pf$nodeA %in% geneInfo_pf$hgnc_symbol,]
  inxDB_pf <- inxDB_pf[inxDB_pf$nodeB %in% geneInfo_pf$hgnc_symbol,]
  
  save(inxDB_pf,geneInfo_pf,file=paste0(netDir,"LRNforCCCnetPrediction_editedbyBI_noECM_noNR.RData"))
  rm(inxDB_pf1)
} else {
  load(paste0(netDir,"LRNforCCCnetPrediction_editedbyBI_noECM_noNR.RData"))
}

# input_LRN_stats ---------------
temp <- mapply(function(A,B) paste(sort(geneInfo_pf[c(A,B),"protein_type"]),collapse="~"),
               A=inxDB_pf$nodeA,B=inxDB_pf$nodeB)
par(mar=c(3,18,1,1),mgp=2:0,las=1)
barplot(sort(table(temp)),horiz=T)
par(mar=c(3,10,1,1),mgp=2:0,las=1)
barplot(sort(table(geneInfo_pf$protein_type)),horiz=T)

# map_to_human --------------------
if (file.exists(paste0(dataPath,".mouse2human.RData"))) {
  load(paste0(dataPath,".mouse2human.RData"))
} else {
  require(biomaRt)
  mart <- useMart("ensembl","mmusculus_gene_ensembl")
  m2h <- getBM(attributes=c("external_gene_name","hsapiens_homolog_associated_gene_name"),
               mart=mart,filters="mgi_symbol",
               values=Reduce(union,sapply(deL,rownames)))
  m2h <- m2h[!m2h$hsapiens_homolog_associated_gene_name == "",]
  save(m2h,file=paste0(dataPath,".mouse2human.RData"))
}
m2h <- m2h[m2h$hsapiens_homolog_associated_gene_name %in% geneInfo_pf$hgnc_symbol,]

## using only uniquely mapping homologs ##
temp <- unique(m2h[duplicated(m2h$hsapiens_homolog_associated_gene_name),
                   "hsapiens_homolog_associated_gene_name"])
m2h <- m2h[!m2h$hsapiens_homolog_associated_gene_name %in% temp,]
temp <- unique(m2h[duplicated(m2h$external_gene_name),"external_gene_name"])
m2h <- m2h[!m2h$external_gene_name %in% temp,]

geneInfo <- geneInfo_pf[geneInfo_pf$hgnc_symbol %in% m2h$hsapiens_homolog_associated_gene_name,]
geneInfo$mgi_symbol <- sapply(geneInfo$hgnc_symbol,function(X)
  m2h$external_gene_name[m2h$hsapiens_homolog_associated_gene_name == X])
rownames(geneInfo) <- geneInfo$mgi_symbol

inxDB <- inxDB_pf[inxDB_pf$nodeA %in% geneInfo$hgnc_symbol & inxDB_pf$nodeB %in% geneInfo$hgnc_symbol,]
inxDB$nodeA <- sapply(inxDB$nodeA,function(X) geneInfo$mgi_symbol[geneInfo$hgnc_symbol == X])
inxDB$nodeB <- sapply(inxDB$nodeB,function(X) geneInfo$mgi_symbol[geneInfo$hgnc_symbol == X])
rownames(inxDB) <- inxDB$key <- paste(inxDB$nodeA,inxDB$nodeB,sep="_")

# data_LRN_stats -------------------
temp <- mapply(function(A,B) paste(sort(geneInfo[c(A,B),"protein_type"]),collapse="~"),
               A=inxDB$nodeA,B=inxDB$nodeB)
par(mar=c(3,18,1,1),mgp=2:0,las=1)
barplot(sort(table(temp)),horiz=T)
par(mar=c(3,10,1,1),mgp=2:0,las=1)
barplot(sort(table(geneInfo$protein_type)),horiz=T)

# ligand_receptor_network -------------------------
edgeThreshold <- 0
sumAbsScoreThreshold <- -log10(0.01)
  
if (!file.exists(paste0(dataPath,".inx_",sumAbsScoreThreshold,"cutoff.RData"))) {
  inx <- inxNode <- list()
  for (a in names(deL)[-length(deL)]) {
    for (b in names(deL)[-(1:which(names(deL) == a))]) {
      l <- paste(a,b,sep="~")
      temp_a <- rownames(deL[[a]])[abs(deL[[a]]$scoreYXoverOX) >= edgeThreshold]
      temp_b <- rownames(deL[[b]])[abs(deL[[b]]$scoreYXoverOX) >= edgeThreshold]
      if (length(temp_a) < 1 | length(temp_b) < 1) { next }
      inx[[l]] <- inxNode[[l]] <- c()
      
      keysAB <- inxDB$key[inxDB$nodeA %in% temp_a & inxDB$nodeB %in% temp_b]
      sumAbsAB <- sapply(strsplit(keysAB,"_"),function(X)
        sum(abs(c(deL[[a]][X[1],"scoreYXoverOX"],deL[[b]][X[2],"scoreYXoverOX"]))))
      sumScaledAB <- sapply(strsplit(keysAB,"_"),function(X)
        sum(c(deL[[a]][X[1],"scoreYXoverOXscaled"],deL[[b]][X[2],"scoreYXoverOXscaled"])))
      if (any(sumAbsAB >= sumAbsScoreThreshold)) {
        keysAB <- keysAB[sumAbsAB >= sumAbsScoreThreshold]
        namesAB <- sapply(strsplit(keysAB,"_"),function(X)
          paste(paste(X[1],a,sep="_"),paste(X[2],b,sep="_"),sep="~"))
        temp_inx <- inxDB[keysAB,c("key","nodeA","nodeB")]
        rownames(temp_inx) <- temp_inx$key <- namesAB
        colnames(temp_inx)[2:3] <- c("geneA","geneB")
        temp_inx$nodeA <- paste(temp_inx$geneA,a,sep="_")
        temp_inx$nodeB <- paste(temp_inx$geneB,b,sep="_")
        temp_inx <- temp_inx[,c(1,4,5,2,3)]
        temp_inx$cellTypeA <- a
        temp_inx$cellTypeB <- b
        temp_inx$proteinTypeA <- geneInfo[temp_inx$geneA,"protein_type"]
        temp_inx$proteinTypeB <- geneInfo[temp_inx$geneB,"protein_type"]
        temp_inx$sumAbsScore <- sumAbsAB[sumAbsAB >= sumAbsScoreThreshold]
        temp_inx$sumScaledScore <- sumScaledAB[sumAbsAB >= sumAbsScoreThreshold]
        inx[[l]] <- rbind(inx[[l]],temp_inx)
      }
      
      keysBA <- inxDB$key[inxDB$nodeA %in% temp_b & inxDB$nodeB %in% temp_a]
      sumAbsBA <- sapply(strsplit(keysBA,"_"),function(X)
        sum(abs(c(deL[[b]][X[1],"scoreYXoverOX"],deL[[a]][X[2],"scoreYXoverOX"]))))
      sumScaledBA <- sapply(strsplit(keysBA,"_"),function(X)
        sum(c(deL[[b]][X[1],"scoreYXoverOXscaled"],deL[[a]][X[2],"scoreYXoverOXscaled"])))
      if (any(sumAbsBA >= sumAbsScoreThreshold)) {
        keysBA <- keysBA[sumAbsBA >= sumAbsScoreThreshold]
        namesBA <- sapply(strsplit(keysBA,"_"),function(X)
          paste(paste(X[2],a,sep="_"),paste(X[1],b,sep="_"),sep="~"))
        temp_inx <- inxDB[keysBA,c("key","nodeB","nodeA")]
        rownames(temp_inx) <- temp_inx$key <- namesBA
        colnames(temp_inx)[2:3] <- c("geneA","geneB")
        temp_inx$nodeA <- paste(temp_inx$geneA,a,sep="_")
        temp_inx$nodeB <- paste(temp_inx$geneB,b,sep="_")
        temp_inx <- temp_inx[,c(1,4,5,2,3)]
        temp_inx$cellTypeA <- a
        temp_inx$cellTypeB <- b
        temp_inx$proteinTypeA <- geneInfo[temp_inx$geneA,"protein_type"]
        temp_inx$proteinTypeB <- geneInfo[temp_inx$geneB,"protein_type"]
        temp_inx$sumAbsScore <- sumAbsBA[sumAbsBA >= sumAbsScoreThreshold]
        temp_inx$sumScaledScore <- sumScaledBA[sumAbsBA >= sumAbsScoreThreshold]
        inx[[l]] <- rbind(inx[[l]],temp_inx)
      }
      if (is.null(inx[[l]])) { next }
      inx[[l]]$sumAbsFactor <- cut(inx[[l]]$sumAbsScore,10,labels=F)
      inx[[l]]$sumScaledFactor <- 6
      if (any(inx[[l]]$sumScaledScore < 0)) {
        inx[[l]]$sumScaledFactor[inx[[l]]$sumScaledScore < 0] <- 
          cut(inx[[l]]$sumScaledScore[inx[[l]]$sumScaledScore < 0],5,labels=F)
      }
      if (any(inx[[l]]$sumScaledScore > 0)) {
        inx[[l]]$sumScaledFactor[inx[[l]]$sumScaledScore > 0] <- 
          inx[[l]]$sumScaledFactor[inx[[l]]$sumScaledScore > 0] +
          cut(inx[[l]]$sumScaledScore[inx[[l]]$sumScaledScore > 0],5,labels=F)
      }
      inx[[l]]$direction <- apply(apply(inx[[l]][,c("proteinTypeA","proteinTypeB")],1,function(X)
        c(grepl("Ligand",X),grepl("Receptor",X))),2,function(Y) {
          if (all(1:4 %in% which(Y))) {
            return("bothLR")
          } else if (all(c(1,4) %in% which(Y))) { # ligand from EC, with receptor in other cell type
            return("LtoR")
          } else if (all(c(2,3) %in% which(Y))) { # ligand from other cell type, receptor in EC
            return("RtoL")
          } else { return("tryECM") } # both ligands or both receptors
        })
      if (any(inx[[l]]$direction == "tryECM")) {
        inx[[l]]$direction[inx[[l]]$direction == "tryECM"] <- apply(apply(inx[[l]][
          inx[[l]]$direction == "tryECM",c("proteinTypeA","proteinTypeB")],1,function(X)
            c(grepl("ECM",X),grepl("Receptor",X))),2,function(Y) {
              if (all(1:4 %in% which(Y))) {
                return("bothER")
              } else if (all(c(1,4) %in% which(Y))) { # ligand from EC, with receptor in other cell type
                return("EtoR")
              } else if (all(c(2,3) %in% which(Y))) { # ligand from other cell type, receptor in EC
                return("RtoE")
              } else { return("tryEL") } # both ligands or both receptors
            })
      }
      if (any(inx[[l]]$direction == "tryEL")) {
        inx[[l]]$direction[inx[[l]]$direction == "tryEL"] <- apply(apply(inx[[l]][
          inx[[l]]$direction == "tryEL",c("proteinTypeA","proteinTypeB")],1,function(X)
            c(grepl("ECM",X),grepl("Ligand",X))),2,function(Y) {
              if (all(1:4 %in% which(Y))) {
                return("bothEL")
              } else if (all(c(1,4) %in% which(Y))) { # ligand from EC, with receptor in other cell type
                return("EtoL")
              } else if (all(c(2,3) %in% which(Y))) { # ligand from other cell type, receptor in EC
                return("LtoE")
              } else { return("None") } # both ligands or both receptors
            })
      }
      
      inxNode[[l]] <- data.frame(node=c(unique(inx[[l]]$nodeA),unique(inx[[l]]$nodeB)),
                                 gene=c(unique(inx[[l]]$geneA),unique(inx[[l]]$geneB)),stringsAsFactors=F)
      rownames(inxNode[[l]]) <- inxNode[[l]]$node
      inxNode[[l]]$cellType <- sapply(strsplit(inxNode[[l]]$node,"_"),function(X) X[2])
      inxNode[[l]]$proteinType <- geneInfo[inxNode[[l]]$gene,"protein_type"]
      inxNode[[l]] <- cbind(inxNode[[l]],
                            rbind(deL[[a]][unique(inx[[l]]$geneA),],
                                  deL[[b]][unique(inx[[l]]$geneB),]))
    }
  }
  save(inx,inxNode,file=paste0(dataPath,".inx_",sumAbsScoreThreshold,"cutoff.RData"))
} else {
  load(paste0(dataPath,".inx_",sumAbsScoreThreshold,"cutoff.RData"))
}

# write_network_tables ----------------
for (l in names(inx)) {
  write.table(inx[[l]][order(inx[[l]]$sumAbsScore,decreasing=T),],
              file=paste0(dataPath,sub("~","_",l),"_",sumAbsScoreThreshold,"cutoff_edges.txt"),
              sep="\t",quote=F,row.names=F,col.names=T)
  write.table(inxNode[[l]][order(abs(inxNode[[l]]$scoreYXoverOX),decreasing=T),],
              file=paste0(dataPath,sub("~","_",l),"_",sumAbsScoreThreshold,"cutoff_nodes.txt"),
              sep="\t",quote=F,row.names=F,col.names=T)
}
temp_inx <- list()
for (l in names(inx)[!grepl("pOPC|cOPC",names(inx))]) {
  if (all(!inx[[l]]$direction != "None")) { next }
  temp_inx[[l]] <- inx[[l]][inx[[l]]$direction != "None",
                       c("key","nodeA","nodeB","geneA","geneB","cellTypeA","cellTypeB",
                         "proteinTypeA","proteinTypeB","sumScaledScore")]
  temp_nodes <- mapply(function(A,B) c(inxNode[[l]][A,"scoreYXoverOX"],
                                       inxNode[[l]][B,"scoreYXoverOX"]),
                       A=temp_inx[[l]]$nodeA,B=temp_inx[[l]]$nodeB)
  temp_inx[[l]]$DEscoreYoungOverOldNodeA <- temp_nodes[1,]
  temp_inx[[l]]$DEscoreYoungOverOldNodeB <- temp_nodes[2,]
  temp_inx[[l]]$scaledSumDEscore <- temp_inx[[l]]$sumScaledScore
  temp_inx[[l]] <- temp_inx[[l]][,-which(colnames(temp_inx[[l]]) == "sumScaledScore")]
}
temp_inxMerge <- do.call(rbind,temp_inx)
write.table(temp_inxMerge,file=paste0(dataPath,".edgeList_",
                                      sumAbsScoreThreshold,"cutoff_ForPub.txt"),
            sep="\t",quote=F,row.names=F,col.names=T)



# plot_network ----------------

prepBipartiteGraphSpread <- function(l,tf) {
  if (!tf %in% c("LtoR","RtoL","RtoE","EtoR","LtoE","EtoL")) {
    stop(paste("tf must be one of:",
               paste(c("LtoR","RtoL","RtoE","EtoR","LtoE","EtoL"),collapse=";"))) 
  }
  co <- strsplit(l,"~")[[1]]
  names(co) <- strsplit(tf,"to")[[1]]

  temp_both <- unique(grep("both",inx[[l]]$direction,value=T))
  temp_both <- temp_both[grepl(names(co)[1],temp_both) & grepl(names(co)[2],temp_both)]
  edl <- inx[[l]]$direction %in% c(tf,temp_both)
  if (!any(edl)) { return(NA) }
  
  temp <- list()
  temp[[co[1]]] <- inxNode[[l]][inx[[l]]$nodeA[edl],
                                c("gene","scoreYXoverOX","scoreYXoverOXscaled","scaledFactor")]
  temp[[co[1]]] <- temp[[co[1]]][!duplicated(temp[[co[1]]]),]
  temp[[co[1]]]$ordered <- rank(
    temp[[co[1]]]$scoreYXoverOX,ties.method="random")/(nrow(temp[[co[1]]]) + 1)*100
  
  temp[[co[2]]] <- inxNode[[l]][inx[[l]]$nodeB[edl],
                                c("gene","scoreYXoverOX","scoreYXoverOXscaled","scaledFactor")]
  temp[[co[2]]] <- temp[[co[2]]][!duplicated(temp[[co[2]]]),]
  temp[[co[2]]]$ordered <- rank(
    temp[[co[2]]]$scoreYXoverOX,ties.method="random")/(nrow(temp[[co[2]]]) + 1)*100

  temp$e <- inx[[l]][edl,c("nodeA","nodeB","sumAbsFactor","sumScaledFactor")]
  colnames(temp$e)[1:2] <- co
  return(temp)
}

plotBipartiteGraphSpread <- function(l,temp) {
  co <- strsplit(l,"~")[[1]]
  names(co) <- strsplit(tf,"to")[[1]]
  co <- co[order(names(co))]
  names(co) <- sapply(names(co),function(X) switch(X,L="Ligands",R="Receptors",E="ECM"))
  temp_main <- paste(names(co)[1],"from",co[1],"to",names(co)[2],"from",co[2])
  
  plot(x=NULL,y=NULL,xlim=c(.1,8.9),ylim=c(0,100),yaxs="i",xaxs="i",
       xaxt="n",yaxt="n",xlab=NA,ylab=NA,main=temp_main,bty="n")
  points(c(rep(2,nrow(temp[[co[1]]])),rep(7,nrow(temp[[co[2]]]))),
         c(temp[[co[1]]]$ordered,temp[[co[2]]]$ordered),pch=21,cex=2,
         col=brewer.pal(11,"PRGn")[c(temp[[co[1]]]$scaledFactor,temp[[co[2]]]$scaledFactor)],
         bg=alpha(brewer.pal(11,"PRGn")[c(temp[[co[1]]]$scaledFactor,temp[[co[2]]]$scaledFactor)],.5))
  
  temp_junk <- apply(temp$e,1,function(X) 
    lines(x=c(2,7),y=c(temp[[co[1]]][X[[co[1]]],"ordered"],
                       temp[[co[2]]][X[[co[2]]],"ordered"]),
          lwd=seq(1,4,length.out=10)[as.integer(X["sumAbsFactor"])],
          col=alpha(brewer.pal(11,"PRGn")[as.integer(X["sumScaledFactor"])],
                    seq(.6,1,length.out=10)[as.integer(X["sumAbsFactor"])])))
  
  text(x=2,y=temp[[co[1]]]$ordered,labels=temp[[co[1]]]$gene,pos=2,col="black")
  text(x=7,y=temp[[co[2]]]$ordered,labels=temp[[co[2]]]$gene,pos=4,col="black")
  
  mtext(text=co,side=1,at=c(2,7),line=.5,font=2)
  mtext(text=c("Older","Differential Gene Expression","Younger"),
        side=2,at=c(10,50,90),line=1.2,font=2,
        col=c(brewer.pal(11,"PRGn")[1],"black",brewer.pal(11,"PRGn")[11]))
  barplot(rep(-.1,100),col=brewer.pal(11,"PRGn")[cut(1:100,11,labels=F)],
          space=0,border=NA,horiz=T,xaxt="n",add=T)
  barplot(rep(.1,100),col=brewer.pal(11,"PRGn")[cut(1:100,11,labels=F)],
          space=0,border=NA,horiz=T,xaxt="n",add=T)
}

for (l in names(inx)) {
  tempL <- sapply(c("LtoR","RtoL","RtoE","EtoR","LtoE","EtoL"),
                  function(tf) prepBipartiteGraphSpread(l,tf),simplify=F)
  tempL <- tempL[sapply(tempL,typeof) == "list"]
  if (length(tempL) == 0) { next }
  temp_figH <- 0.20 * max(sapply(tempL,function(X) sapply(X[1:2],nrow)))
  if (temp_figH < 5) { temp_figH <- 5 } 
  for (tf in names(tempL)) {
    pdf(file=paste0(dataPath,sub("~","_",l),"_",tf,"_",sumAbsScoreThreshold,"cutoff.pdf"),
        width=6,height=temp_figH)
    par(mar=c(2,3,2,1),mgp=c(1,1,0))
    plotBipartiteGraphSpread(l,tempL[[tf]])
    dev.off()
  }
}


# SumFig_prep --------
CTOI <- c("OPC","OLG","OEG","ASC","mNEUR","EPC","EC","PC","ABC","MG","MAC")
inx <- inx[sapply(strsplit(names(inx),"~"),function(X) all(X %in% CTOI))]
inxNode <- inxNode[sapply(strsplit(names(inxNode),"~"),function(X) all(X %in% CTOI))]
LtoR <- sapply(inx,function(X) sum(X$direction == "LtoR"))
RtoL <- sapply(inx,function(X) sum(X$direction == "RtoL"))


```{r SumFig_js,include=F}
library(networkD3)
links <- do.call(rbind,lapply(strsplit(names(LtoR),"~"),function(X) 
  c(X[1],X[2],LtoR[paste(X,collapse="~")],use.names=F)))
colnames(links) <- c("source","target","value")
links2 <- do.call(rbind,lapply(strsplit(names(RtoL),"~"),function(X) 
  c(X[2],X[1],RtoL[paste(X,collapse="~")],use.names=F)))
colnames(links2) <- c("source","target","value")
links <- as.data.frame(rbind(links,links2))
rm(links2)

nodes <- data.frame(name=rep(CTOI,2))
links$IDsource <- match(links$source,CTOI) - 1 
links$IDtarget <- match(links$target,CTOI) + length(CTOI) - 1


sankeyNetwork(Links = links, Nodes = nodes,
              Source = "IDsource", Target = "IDtarget",
              Value = "value", NodeID = "name", 
              sinksRight=T,iterations=100)
```


```{r SumFig_Rplot,include=F}
library(riverplot)
links <- do.call(rbind,lapply(strsplit(names(LtoR),"~"),function(X) 
  c(X[1],X[2],LtoR[paste(X,collapse="~")],use.names=F)))
links2 <- do.call(rbind,lapply(strsplit(names(RtoL),"~"),function(X) 
  c(X[2],X[1],RtoL[paste(X,collapse="~")],use.names=F)))
links <- as.data.frame(rbind(links,links2),stringsAsFactors=F)
rm(links2)
colnames(links) <- c("N1","N2","Value")
links$N1 <- paste0(links$N1,"_ligand")
links$N2 <- paste0(links$N2,"_receptor")
links$Value <- as.integer(links$Value)
links <- links[links$Value > 0,]
links$lty <- 0

nodes <- data.frame(ID=c(paste0(CTOI,"_ligand"),paste0(CTOI,"_receptor")),
                    x=c(rep(1,length(CTOI)),rep(2,length(CTOI))),
                    labels=rep(CTOI,2),
                    col=rep(c("olivedrab4","olivedrab3","olivedrab1","steelblue1",
                              "magenta3","lightgoldenrod4","sienna4","sienna3",
                              "peachpuff3","red3","red1"),2),
                    stringsAsFactors=F)
rownames(nodes) <- nodes$ID
nodes <- nodes[nodes$ID %in% unique(unlist(links[,1:2])),]
nodes$srt <- 0
nodes$lty <- 1
nodes <- nodes[c(paste0(c("PC","MG","EC","OLG","ASC","mNEUR",
                          "MAC","EPC","OEG","OPC","ABC"),"_ligand"),
                 paste0(c("MAC","EPC","mNEUR","ASC","OEG","ABC",
                          "OPC","PC","MG","EC","OLG"),"_receptor")),]

SD <- makeRiver(nodes,links)

pdf(file=paste0(dataPath,"summary_",sumAbsScoreThreshold,"cutoff.pdf"),
    width=6,height=8)
par(mar=rep(1,4),mgp=2:0)
riverplot(SD,plot_area=1,fix.pdf=T)
mtext("Ligands",side=2,font=2,cex=1.5,line=-.5)
mtext("Receptors",side=4,font=2,cex=1.5,line=-.5)
dev.off()

```


