# library(stringr)
# library("protr")
# library("bio3d")
# library(dplyr)

######## Normalized table #######
normalizeDF <- function(x){
  # x <- na.omit(x)
  rs <- x
  for(i in 1:ncol(x)){
    tt <- normalize(x[,i])
    rs[,i] <- tt[1,]
  }
  rs
}

########################## Amino Acid Composition (AAC) ##########################
generateAAC <- function(dt){
  seq_name <- c()
  seq_length <- c()
  label <- c()
  rs <- data.frame()
  for(n in 1:nrow(dt)){
    aac <- extractAAC(dt[n, "seq_name"])
    df <- as.data.frame(split(unname(aac), sub(".*[.]", "", names(aac))), stringsAsFactors = FALSE)
    # df$label <- dt[n, "label"]
    rs <- rbind(rs, df)
  }
  rs
}

########################## Pseudo Amino Acid Composition (PAAC) ##########################
generatePAAC <- function(dt){
  seq_name <- c()
  seq_length <- c()
  label <- c()
  rs <- data.frame()
  for(n in 1:nrow(dt)){
    paac <- extractPAAC(dt[n, "seq_name"], lambda = 4)
    df <- as.data.frame(split(unname(paac), sub(".*[.]", "", names(paac))), stringsAsFactors = FALSE)
    # df$label <- dt[n, "label"]
    rs <- rbind(rs, df)
  }
  rs
}

########################## Quasi Sequence Order Descriptor (QSO) ########################## 
generateQSO <- function(dt){
  seq_name <- c()
  seq_length <- c()
  label <- c()
  rs <- data.frame()
  for(n in 1:nrow(dt)){
    qso <- extractQSO(dt[n, "seq_name"], nlag = 4)
    df <- as.data.frame(split(unname(qso), sub("*[.]", ".", names(qso))), stringsAsFactors = FALSE)
    # df$label <- dt[n, "label"]
    rs <- rbind(rs, df)
  }
  rs
}

######################################### AAindex ######################################### 
# impIndexGroup <- c("HUTJ700101","SUEM840102","GEIM800104","CHAM830105","ARGP820103", "QIAN880128","ARGP820102","COSI940101","PALJ810109")
df.aaindex <- read.csv('aaindex.csv') # 552 indices
df.aaindex <- df.aaindex[,2:11]
##################### AAindex with 7 types of group ######################
hydrophobicity.g1 <- c('R', 'K', 'E', 'D', 'Q', 'N')
hydrophobicity.g2 <- c('G', 'A', 'S', 'T', 'P', 'H', 'Y')
hydrophobicity.g3 <- c('C', 'L', 'V', 'I', 'M', 'F', 'W')

vanderwa.g1 <- c('G', 'A', 'S', 'T', 'P', 'D', 'C')
vanderwa.g2 <- c('N', 'V', 'E', 'Q', 'I', 'L')
vanderwa.g3 <- c('M', 'H', 'K', 'F', 'R', 'Y', 'W')

polarity.g1 <- c('L', 'I', 'F', 'W', 'C', 'M', 'V', 'Y')
polarity.g2 <- c('P', 'A', 'T', 'G', 'S')
polarity.g3 <- c('H', 'Q', 'R', 'K', 'N', 'E', 'D')

polarizability.g1 <- c('G', 'A', 'S', 'D', 'T')
polarizability.g2 <- c('C', 'P', 'N', 'V', 'E', 'Q', 'I', 'L')
polarizability.g3 <- c('K', 'M', 'H', 'F', 'R', 'Y', 'W')

charge.g1 <- c('K', 'R')
charge.g2 <- c('A', 'N', 'C', 'Q', 'G', 'H', 'I', 'L', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')
charge.g3 <- c('D', 'E')

second.struct.g1 <- c('E', 'A', 'L', 'M', 'Q', 'K', 'R', 'H')
second.struct.g2 <- c('V', 'I', 'Y', 'C', 'W', 'F', 'T')
second.struct.g3 <- c('G', 'N', 'P', 'S', 'D')

solv.acc.g1 <- c('A', 'L', 'F', 'C', 'G', 'I', 'V', 'W')
solv.acc.g2 <- c('R', 'K', 'Q', 'E', 'N', 'D')
solv.acc.g3 <- c('M', 'S', 'P', 'T', 'H', 'Y')

################ can paste multiple indices #################################
################ calculate 3 properties: hydro,vande, charg #################
generateAAindexWithWeight_NF <- function(df, ind){
  df.rs <- data.frame(seq_name = df$seq_name)
  df.aaindex.local <- df.aaindex[, c('aa', ind)]
  for(i in 2:ncol(df.aaindex.local)){
    ######## calculate each properties of each index ########
    df.sub <- generateAverageWeight_NF(df,colnames(df.aaindex.local)[i])
    df.rs <- cbind(df.rs, df.sub)
  }
  df.rs <- df.rs %>% select(-seq_name)
  df.rs
}
################ calculate 4 properties: polty,polzt,secst,solva #################
generateAAindexWithWeight_NF2 <- function(df, ind){
  df.rs <- data.frame(seq_name = df$seq_name)
  df.aaindex.local <- df.aaindex[, c('aa', ind)]
  for(i in 2:ncol(df.aaindex.local)){
    df.sub <- generateAverageWeight_NF2(df,colnames(df.aaindex.local)[i])
    df.rs <- cbind(df.rs, df.sub)
  }
  df.rs <- df.rs %>% select(-seq_name)
  df.rs
}

################ calculate 3 properties: hydro,vande, charg #################
generateAverageWeight_NF <- function(df, ind){
  hydro.g1 <- c()
  hydro.g2 <- c()
  hydro.g3 <- c()
  vande.g1 <- c()
  vande.g2 <- c()
  vande.g3 <- c()
  charg.g1 <- c()
  charg.g2 <- c()
  charg.g3 <- c()
  
  for(i in 1:nrow(df)){
    hydro.g1.freq <- data.frame(aa = hydrophobicity.g1, qty=str_count(df[i,]$seq_name, hydrophobicity.g1))#count each member of each group ex: group1 of hydro
    hydro.g1.aaindex <- df.aaindex[df.aaindex$aa %in% hydrophobicity.g1, c('aa', ind)]# AAindex values of member group 1
    hydro.g1.aaindex <- merge(hydro.g1.freq, hydro.g1.aaindex, by = 'aa') ## aa, count, index_val of group 1
    ######### Sum frequencies x average_weight (exclude non-exist member) / seq_length
    hydro.g1 <- c(hydro.g1, sum(hydro.g1.aaindex[hydro.g1.aaindex$qty>0,]$qty)*mean(hydro.g1.aaindex[hydro.g1.aaindex$qty>0,c(ind)])/df[i,]$seq_length)
    hydro.g2.freq <- data.frame(aa = hydrophobicity.g2, qty=str_count(df[i,]$seq_name, hydrophobicity.g2))
    hydro.g2.aaindex <- df.aaindex[df.aaindex$aa %in% hydrophobicity.g2, c('aa', ind)]
    hydro.g2.aaindex <- merge(hydro.g2.freq, hydro.g2.aaindex, by = 'aa')
    hydro.g2 <- c(hydro.g2, sum(hydro.g2.aaindex[hydro.g2.aaindex$qty>0,]$qty)*mean(hydro.g2.aaindex[hydro.g2.aaindex$qty>0,c(ind)])/df[i,]$seq_length)
    hydro.g3.freq <- data.frame(aa = hydrophobicity.g3, qty=str_count(df[i,]$seq_name, hydrophobicity.g3))
    hydro.g3.aaindex <- df.aaindex[df.aaindex$aa %in% hydrophobicity.g3, c('aa', ind)]
    hydro.g3.aaindex <- merge(hydro.g3.freq, hydro.g3.aaindex, by = 'aa')
    hydro.g3 <- c(hydro.g3, sum(hydro.g3.aaindex[hydro.g3.aaindex$qty>0,]$qty)*mean(hydro.g3.aaindex[hydro.g3.aaindex$qty>0,c(ind)])/df[i,]$seq_length)
    
    vande.g1.freq <- data.frame(aa = vanderwa.g1, qty = str_count(df[i,]$seq_name,vanderwa.g1))
    vande.g1.aaindex <- df.aaindex[df.aaindex$aa %in%vanderwa.g1, c('aa', ind)]
    vande.g1.aaindex <- merge(vande.g1.freq, vande.g1.aaindex, by = 'aa')
    vande.g1 <- c(vande.g1, sum(vande.g1.aaindex[vande.g1.aaindex$qty>0,]$qty)*mean(vande.g1.aaindex[vande.g1.aaindex$qty>0,c(ind)])/df[i,]$seq_length)
    vande.g2.freq <- data.frame(aa =vanderwa.g2, qty=str_count(df[i,]$seq_name,vanderwa.g2))
    vande.g2.aaindex <- df.aaindex[df.aaindex$aa %in%vanderwa.g2, c('aa', ind)]
    vande.g2.aaindex <- merge(vande.g2.freq, vande.g2.aaindex, by = 'aa')
    vande.g2 <- c(vande.g2, sum(vande.g2.aaindex[vande.g2.aaindex$qty>0,]$qty)*mean(vande.g2.aaindex[vande.g2.aaindex$qty>0,c(ind)])/df[i,]$seq_length)
    vande.g3.freq <- data.frame(aa =vanderwa.g3, qty=str_count(df[i,]$seq_name,vanderwa.g3))
    vande.g3.aaindex <- df.aaindex[df.aaindex$aa %in%vanderwa.g3, c('aa', ind)]
    vande.g3.aaindex <- merge(vande.g3.freq, vande.g3.aaindex, by = 'aa')
    vande.g3 <- c(vande.g3, sum(vande.g3.aaindex[vande.g3.aaindex$qty>0,]$qty)*mean(vande.g3.aaindex[vande.g3.aaindex$qty>0,c(ind)])/df[i,]$seq_length)
    
    charg.g1.freq <- data.frame(aa =charge.g1, qty=str_count(df[i,]$seq_name,charge.g1))
    charg.g1.aaindex <- df.aaindex[df.aaindex$aa %in%charge.g1, c('aa', ind)]
    charg.g1.aaindex <- merge(charg.g1.freq, charg.g1.aaindex, by = 'aa')
    charg.g1 <- c(charg.g1, sum(charg.g1.aaindex[charg.g1.aaindex$qty>0,]$qty)*mean(charg.g1.aaindex[charg.g1.aaindex$qty>0,c(ind)])/df[i,]$seq_length)
    charg.g2.freq <- data.frame(aa =charge.g2, qty=str_count(df[i,]$seq_name,charge.g2))
    charg.g2.aaindex <- df.aaindex[df.aaindex$aa %in%charge.g2, c('aa', ind)]
    charg.g2.aaindex <- merge(charg.g2.freq, charg.g2.aaindex, by = 'aa')
    charg.g2 <- c(charg.g2, sum(charg.g2.aaindex[charg.g2.aaindex$qty>0,]$qty)*mean(charg.g2.aaindex[charg.g2.aaindex$qty>0,c(ind)])/df[i,]$seq_length)
    charg.g3.freq <- data.frame(aa =charge.g3, qty=str_count(df[i,]$seq_name,charge.g3))
    charg.g3.aaindex <- df.aaindex[df.aaindex$aa %in%charge.g3, c('aa', ind)]
    charg.g3.aaindex <- merge(charg.g3.freq, charg.g3.aaindex, by = 'aa')
    charg.g3 <- c(charg.g3, sum(charg.g3.aaindex[charg.g3.aaindex$qty>0,]$qty)*mean(charg.g3.aaindex[charg.g3.aaindex$qty>0,c(ind)])/df[i,]$seq_length)
  }
  
  df.rs <- data.frame(
    hydro.g1 = hydro.g1,
    hydro.g2 = hydro.g2,
    hydro.g3 = hydro.g3,
    vande.g1 = vande.g1,
    vande.g2 = vande.g2,
    vande.g3 = vande.g3,
    charg.g1 = charg.g1,
    charg.g2 = charg.g2,
    charg.g3 = charg.g3
    
  )
  colnames(df.rs) <- c(paste(ind, 'hydrophobicity.nf.g1', sep='.'),paste(ind, 'hydrophobicity.nf.g2', sep='.'),paste(ind, 'hydrophobicity.nf.g3', sep='.')
                       ,paste(ind, 'vanderwa.nf.g1', sep='.'),paste(ind, 'vanderwa.nf.g2', sep='.'),paste(ind, 'vanderwa.nf.g3', sep='.')
                       ,paste(ind, 'charge.nf.g1', sep='.'),paste(ind, 'charge.nf.g2', sep='.'),paste(ind, 'charge.nf.g3', sep='.'))
  
  df.rs[is.na(df.rs)]=0
  #df.rs <- normalizeDF(select(df.rs, -seq_name))
  df.rs
}

################ calculate 4 properties: polty,polzt,secst,solva #################
generateAverageWeight_NF2 <- function(df, ind){
  
  polty.g1 <- c()
  polty.g2 <- c()
  polty.g3 <- c()
  
  polzt.g1 <- c()
  polzt.g2 <- c()
  polzt.g3 <- c()
  
  secst.g1 <- c()
  secst.g2 <- c()
  secst.g3 <- c()
  
  solva.g1 <- c()
  solva.g2 <- c()
  solva.g3 <- c()
  
  for(i in 1:nrow(df)){
    polty.g1.freq <- data.frame(aa =polarity.g1, qty=str_count(df[i,]$seq_name,polarity.g1))
    polty.g1.aaindex <- df.aaindex[df.aaindex$aa %in%polarity.g1, c('aa', ind)]
    polty.g1.aaindex <- merge(polty.g1.freq, polty.g1.aaindex, by = 'aa')
    polty.g1 <- c(polty.g1, sum(polty.g1.aaindex[polty.g1.aaindex$qty>0,]$qty)*mean(polty.g1.aaindex[polty.g1.aaindex$qty>0,c(ind)])/df[i,]$seq_length)
    polty.g2.freq <- data.frame(aa =polarity.g2, qty=str_count(df[i,]$seq_name,polarity.g2))
    polty.g2.aaindex <- df.aaindex[df.aaindex$aa %in%polarity.g2, c('aa', ind)]
    polty.g2.aaindex <- merge(polty.g2.freq, polty.g2.aaindex, by = 'aa')
    polty.g2 <- c(polty.g2, sum(polty.g2.aaindex[polty.g2.aaindex$qty>0,]$qty)*mean(polty.g2.aaindex[polty.g2.aaindex$qty>0,c(ind)])/df[i,]$seq_length)
    polty.g3.freq <- data.frame(aa =polarity.g3, qty=str_count(df[i,]$seq_name,polarity.g3))
    polty.g3.aaindex <- df.aaindex[df.aaindex$aa %in%polarity.g3, c('aa', ind)]
    polty.g3.aaindex <- merge(polty.g3.freq, polty.g3.aaindex, by = 'aa')
    polty.g3 <- c(polty.g3, sum(polty.g3.aaindex[polty.g3.aaindex$qty>0,]$qty)*mean(polty.g3.aaindex[polty.g3.aaindex$qty>0,c(ind)])/df[i,]$seq_length)
    
    polzt.g1.freq <- data.frame(aa =polarizability.g1, qty=str_count(df[i,]$seq_name,polarizability.g1))
    polzt.g1.aaindex <- df.aaindex[df.aaindex$aa %in%polarizability.g1, c('aa', ind)]
    polzt.g1.aaindex <- merge(polzt.g1.freq, polzt.g1.aaindex, by = 'aa')
    polzt.g1 <- c(polzt.g1, sum(polzt.g1.aaindex[polzt.g1.aaindex$qty>0,]$qty)*mean(polzt.g1.aaindex[polzt.g1.aaindex$qty>0,c(ind)])/df[i,]$seq_length)
    polzt.g2.freq <- data.frame(aa =polarizability.g2, qty=str_count(df[i,]$seq_name,polarizability.g2))
    polzt.g2.aaindex <- df.aaindex[df.aaindex$aa %in%polarizability.g2, c('aa', ind)]
    polzt.g2.aaindex <- merge(polzt.g2.freq, polzt.g2.aaindex, by = 'aa')
    polzt.g2 <- c(polzt.g2, sum(polzt.g2.aaindex[polzt.g2.aaindex$qty>0,]$qty)*mean(polzt.g2.aaindex[polzt.g2.aaindex$qty>0,c(ind)])/df[i,]$seq_length)
    polzt.g3.freq <- data.frame(aa =polarizability.g3, qty=str_count(df[i,]$seq_name,polarizability.g3))
    polzt.g3.aaindex <- df.aaindex[df.aaindex$aa %in%polarizability.g3, c('aa', ind)]
    polzt.g3.aaindex <- merge(polzt.g3.freq, polzt.g3.aaindex, by = 'aa')
    polzt.g3 <- c(polzt.g3, sum(polzt.g3.aaindex[polzt.g3.aaindex$qty>0,]$qty)*mean(polzt.g3.aaindex[polzt.g3.aaindex$qty>0,c(ind)])/df[i,]$seq_length)
    
    secst.g1.freq <- data.frame(aa =second.struct.g1, qty=str_count(df[i,]$seq_name,second.struct.g1))
    secst.g1.aaindex <- df.aaindex[df.aaindex$aa %in%second.struct.g1, c('aa', ind)]
    secst.g1.aaindex <- merge(secst.g1.freq, secst.g1.aaindex, by = 'aa')
    secst.g1 <- c(secst.g1, sum(secst.g1.aaindex[secst.g1.aaindex$qty>0,]$qty)*mean(secst.g1.aaindex[secst.g1.aaindex$qty>0,c(ind)])/df[i,]$seq_length)
    secst.g2.freq <- data.frame(aa =second.struct.g2, qty=str_count(df[i,]$seq_name,second.struct.g2))
    secst.g2.aaindex <- df.aaindex[df.aaindex$aa %in%second.struct.g2, c('aa', ind)]
    secst.g2.aaindex <- merge(secst.g2.freq, secst.g2.aaindex, by = 'aa')
    secst.g2 <- c(secst.g2, sum(secst.g2.aaindex[secst.g2.aaindex$qty>0,]$qty)*mean(secst.g2.aaindex[secst.g2.aaindex$qty>0,c(ind)])/df[i,]$seq_length)
    secst.g3.freq <- data.frame(aa =second.struct.g3, qty=str_count(df[i,]$seq_name,second.struct.g3))
    secst.g3.aaindex <- df.aaindex[df.aaindex$aa %in%second.struct.g3, c('aa', ind)]
    secst.g3.aaindex <- merge(secst.g3.freq, secst.g3.aaindex, by = 'aa')
    secst.g3 <- c(secst.g3, sum(secst.g3.aaindex[secst.g3.aaindex$qty>0,]$qty)*mean(secst.g3.aaindex[secst.g3.aaindex$qty>0,c(ind)])/df[i,]$seq_length)
    
    solva.g1.freq <- data.frame(aa =solv.acc.g1, qty=str_count(df[i,]$seq_name,solv.acc.g1))
    solva.g1.aaindex <- df.aaindex[df.aaindex$aa %in%solv.acc.g1, c('aa', ind)]
    solva.g1.aaindex <- merge(solva.g1.freq, solva.g1.aaindex, by = 'aa')
    solva.g1 <- c(solva.g1, sum(solva.g1.aaindex[solva.g1.aaindex$qty>0,]$qty)*mean(solva.g1.aaindex[solva.g1.aaindex$qty>0,c(ind)])/df[i,]$seq_length)
    solva.g2.freq <- data.frame(aa =solv.acc.g2, qty=str_count(df[i,]$seq_name,solv.acc.g2))
    solva.g2.aaindex <- df.aaindex[df.aaindex$aa %in%solv.acc.g2, c('aa', ind)]
    solva.g2.aaindex <- merge(solva.g2.freq, solva.g2.aaindex, by = 'aa')
    solva.g2 <- c(solva.g2, sum(solva.g2.aaindex[solva.g2.aaindex$qty>0,]$qty)*mean(solva.g2.aaindex[solva.g2.aaindex$qty>0,c(ind)])/df[i,]$seq_length)
    solva.g3.freq <- data.frame(aa =solv.acc.g3, qty=str_count(df[i,]$seq_name,solv.acc.g3))
    solva.g3.aaindex <- df.aaindex[df.aaindex$aa %in%solv.acc.g3, c('aa', ind)]
    solva.g3.aaindex <- merge(solva.g3.freq, solva.g3.aaindex, by = 'aa')
    solva.g3 <- c(solva.g3, sum(solva.g3.aaindex[solva.g3.aaindex$qty>0,]$qty)*mean(solva.g3.aaindex[solva.g3.aaindex$qty>0,c(ind)])/df[i,]$seq_length)
    
  }
  df.rs <- data.frame(
    
    polty.g1 = polty.g1,
    polty.g2 = polty.g2,
    polty.g3 = polty.g3,
    polzt.g1 = polzt.g1,
    polzt.g2 = polzt.g2,
    polzt.g3 = polzt.g3,
    
    secst.g1 = secst.g1,
    secst.g2 = secst.g2,
    secst.g3 = secst.g3,
    solva.g1 = solva.g1,
    solva.g2 = solva.g2,
    solva.g3 = solva.g3
    
  )
  colnames(df.rs) <- c(paste(ind, 'polarity.nf.g1', sep='.'),paste(ind, 'polarity.nf.g2', sep='.'),paste(ind, 'polarity.nf.g3', sep='.'),paste(ind, 'polarizability.nf.g1', sep='.'),paste(ind, 'polarizability.nf.g2', sep='.'),paste(ind, 'polarizability.nf.g3', sep='.'),
                       paste(ind, 'second.struct.nf.g1', sep='.'),paste(ind, 'second.struct.nf.g2', sep='.'),paste(ind, 'second.struct.nf.g3', sep='.'),paste(ind, 'solv.acc.nf.g1', sep='.'),paste(ind, 'solv.acc.nf.g2', sep='.'),paste(ind, 'solv.acc.nf.g3', sep='.')
  )
  df.rs[is.na(df.rs)]=0
  #df.rs <- normalizeDF(df.rs)
  df.rs
}

################################################ End AAindex ################################ 

################ Generate original features AAC+PAAC+QSO #################  
generate_features <- function(df){
  df.aac <- generateAAC(df)
  df.paac <- generatePAAC(df)
  colnames(df.paac) <- paste("paac", colnames(df.paac), sep = "_")
  # names(df.paac)[25] <- "label"
  df.qso <- generateQSO(df)
  # df.cb <- cbind(select(df.paac, -label), select(df.qso, -label), df.aac)
  df.cb <- cbind(df.paac, df.qso, df.aac)
  df.cb$label <- df$label
  df.rs <- df.cb
  df.rs
}
encode_features <- function(df){
  df.aac <- generateAAC(df)
  df.paac <- generatePAAC(df)
  colnames(df.paac) <- paste("paac", colnames(df.paac), sep = "_")
  df.qso <- generateQSO(df)
  df.cb <- cbind(df.paac, df.qso, df.aac)
  df.rs <- df.cb
  df.rs
}
################ Generate seperate AAindex properties #################    
####### Generate feature 1 AAIndex 001-100 ##################
generate_features11 <- function(df){
  
  df.rs <- data.frame(seq_name = df$seq_name)
  for(i in 2:51){
    df.index.nf1 <- generateAAindexWithWeight_NF(df, colnames(df.aaindex)[i])
    df.index.nf2 <- generateAAindexWithWeight_NF2(df, colnames(df.aaindex)[i])
    df.rs <- cbind(df.rs, df.index.nf1, df.index.nf2)
  }
  df.rs <- df.rs %>% select(-seq_name)
  df.rs
}
generate_features12 <- function(df){
  
  df.rs <- data.frame(seq_name = df$seq_name)
  for(i in 52:101){
    df.index.nf1 <- generateAAindexWithWeight_NF(df, colnames(df.aaindex)[i])
    df.index.nf2 <- generateAAindexWithWeight_NF2(df, colnames(df.aaindex)[i])
    df.rs <- cbind(df.rs, df.index.nf1, df.index.nf2)
  }
  df.rs <- df.rs %>% select(-seq_name)
  df.rs
}
####### Generate feature 2 AAIndex 101-200 ##################
generate_features21 <- function(df){
  
  df.rs <- data.frame(seq_name = df$seq_name)
  for(i in 102:151){
    df.index.nf1 <- generateAAindexWithWeight_NF(df, colnames(df.aaindex)[i])
    df.index.nf2 <- generateAAindexWithWeight_NF2(df, colnames(df.aaindex)[i])
    df.rs <- cbind(df.rs, df.index.nf1, df.index.nf2)
  }
  df.rs <- df.rs %>% select(-seq_name)
  df.rs
}
generate_features22 <- function(df){
  
  df.rs <- data.frame(seq_name = df$seq_name)
  for(i in 152:201){
    df.index.nf1 <- generateAAindexWithWeight_NF(df, colnames(df.aaindex)[i])
    df.index.nf2 <- generateAAindexWithWeight_NF2(df, colnames(df.aaindex)[i])
    df.rs <- cbind(df.rs, df.index.nf1, df.index.nf2)
  }
  df.rs <- df.rs %>% select(-seq_name)
  df.rs
}
####### Generate feature 3 AAIndex 201-300 ##################
generate_features31 <- function(df){
  
  df.rs <- data.frame(seq_name = df$seq_name)
  for(i in 202:251){
    df.index.nf1 <- generateAAindexWithWeight_NF(df, colnames(df.aaindex)[i])
    df.index.nf2 <- generateAAindexWithWeight_NF2(df, colnames(df.aaindex)[i])
    df.rs <- cbind(df.rs, df.index.nf1, df.index.nf2)
  }
  df.rs <- df.rs %>% select(-seq_name)
  df.rs
}
generate_features32 <- function(df){
  
  df.rs <- data.frame(seq_name = df$seq_name)
  for(i in 252:301){
    df.index.nf1 <- generateAAindexWithWeight_NF(df, colnames(df.aaindex)[i])
    df.index.nf2 <- generateAAindexWithWeight_NF2(df, colnames(df.aaindex)[i])
    df.rs <- cbind(df.rs, df.index.nf1, df.index.nf2)
  }
  df.rs <- df.rs %>% select(-seq_name)
  df.rs
}
####### Generate feature 4 AAIndex 301-400 ##################
generate_features41 <- function(df){
  
  df.rs <- data.frame(seq_name = df$seq_name)
  for(i in 302:351){
    df.index.nf1 <- generateAAindexWithWeight_NF(df, colnames(df.aaindex)[i])
    df.index.nf2 <- generateAAindexWithWeight_NF2(df, colnames(df.aaindex)[i])
    df.rs <- cbind(df.rs, df.index.nf1, df.index.nf2)
  }
  df.rs <- df.rs %>% select(-seq_name)
  df.rs
}
generate_features42 <- function(df){
  
  df.rs <- data.frame(seq_name = df$seq_name)
  for(i in 352:401){
    df.index.nf1 <- generateAAindexWithWeight_NF(df, colnames(df.aaindex)[i])
    df.index.nf2 <- generateAAindexWithWeight_NF2(df, colnames(df.aaindex)[i])
    df.rs <- cbind(df.rs, df.index.nf1, df.index.nf2)
  }
  df.rs <- df.rs %>% select(-seq_name)
  df.rs
}
####### Generate feature 5 AAIndex 401-500 ##################
generate_features51 <- function(df){
  
  df.rs <- data.frame(seq_name = df$seq_name)
  for(i in 402:451){
    df.index.nf1 <- generateAAindexWithWeight_NF(df, colnames(df.aaindex)[i])
    df.index.nf2 <- generateAAindexWithWeight_NF2(df, colnames(df.aaindex)[i])
    df.rs <- cbind(df.rs, df.index.nf1, df.index.nf2)
  }
  df.rs <- df.rs %>% select(-seq_name)
  df.rs
}
generate_features52 <- function(df){
  
  df.rs <- data.frame(seq_name = df$seq_name)
  for(i in 452:501){
    df.index.nf1 <- generateAAindexWithWeight_NF(df, colnames(df.aaindex)[i])
    df.index.nf2 <- generateAAindexWithWeight_NF2(df, colnames(df.aaindex)[i])
    df.rs <- cbind(df.rs, df.index.nf1, df.index.nf2)
  }
  df.rs <- df.rs %>% select(-seq_name)
  df.rs
}
####### Generate feature 6 AAIndex 501-553 ##################
generate_features6 <- function(df){
  
  df.rs <- data.frame(seq_name = df$seq_name)
  for(i in 502:553){
    df.index.nf1 <- generateAAindexWithWeight_NF(df, colnames(df.aaindex)[i])
    df.index.nf2 <- generateAAindexWithWeight_NF2(df, colnames(df.aaindex)[i])
    df.rs <- cbind(df.rs, df.index.nf1, df.index.nf2)
  }
  df.rs <- df.rs %>% select(-seq_name)
  df.rs
}
####### Generate only important index properties ##################
generate_impFeaturesOld <- function(df, importantFeatures){
  df.rs <- data.frame(seq_name = df$seq_name)
  for(i in 1:length(importantFeatures)){
    df.index.nf1 <- generateAAindexWithWeight_NF(df, importantFeatures[i])
    df.index.nf2 <- generateAAindexWithWeight_NF2(df, importantFeatures[i])
    df.rs <- cbind(df.rs, df.index.nf1, df.index.nf2)
  }
  df.rs <- df.rs %>% select(-seq_name)
  df.rs
}
#####################
generate_feature1 <- function(df.input){
  df.output1 <- generate_features11(df.input)
  df.output2 <- generate_features12(df.input)
  df.output <- cbind(df.output1,df.output2) 
  df.output
}
generate_feature23 <- function(df.input){
  df.output1 <- generate_features21(df.input)
  df.output2 <- generate_features22(df.input)
  df.output3 <- generate_features31(df.input)
  df.output4 <- generate_features32(df.input)
  df.output <- cbind(df.output1,df.output2,df.output3,df.output4) 
  df.output
}
generate_feature456 <- function(df.input){
  df.output1 <- generate_features41(df.input)
  df.output2 <- generate_features42(df.input)
  df.output3 <- generate_features51(df.input)
  df.output4 <- generate_features52(df.input)
  df.output5 <- generate_features6(df.input)
  df.output <- cbind(df.output1,df.output2,df.output3,df.output4,df.output5) 
  df.output
}

######### Calculate MCC ########
calculateMCC <- function(conf.table){
  MCC <- (conf.table$table[1,1]*conf.table$table[2,2]-conf.table$table[1,2]*conf.table$table[2,1])/(sqrt(sum(conf.table$table[2,])*sum(conf.table$table[,2]))*sqrt(sum(conf.table$table[1,])*sum(conf.table$table[,1])))
  MCC
}

###################### New update calculate 7 properties: hydro,vande, charg #################
output.columns <- c("HUTJ700101.hydrophobicity.nf.g1","HUTJ700101.hydrophobicity.nf.g2","HUTJ700101.hydrophobicity.nf.g3","HUTJ700101.vanderwa.nf.g1","HUTJ700101.vanderwa.nf.g2","HUTJ700101.vanderwa.nf.g3","HUTJ700101.charge.nf.g1","HUTJ700101.charge.nf.g2","HUTJ700101.charge.nf.g3","HUTJ700101.polarity.nf.g1","HUTJ700101.polarity.nf.g2","HUTJ700101.polarity.nf.g3","HUTJ700101.polarizability.nf.g1","HUTJ700101.polarizability.nf.g2","HUTJ700101.polarizability.nf.g3","HUTJ700101.second.struct.nf.g1","HUTJ700101.second.struct.nf.g2","HUTJ700101.second.struct.nf.g3","HUTJ700101.solv.acc.nf.g1","HUTJ700101.solv.acc.nf.g2","HUTJ700101.solv.acc.nf.g3","SUEM840102.hydrophobicity.nf.g1","SUEM840102.hydrophobicity.nf.g2","SUEM840102.hydrophobicity.nf.g3","SUEM840102.vanderwa.nf.g1","SUEM840102.vanderwa.nf.g2","SUEM840102.vanderwa.nf.g3","SUEM840102.charge.nf.g1","SUEM840102.charge.nf.g2","SUEM840102.charge.nf.g3","SUEM840102.polarity.nf.g1","SUEM840102.polarity.nf.g2","SUEM840102.polarity.nf.g3","SUEM840102.polarizability.nf.g1","SUEM840102.polarizability.nf.g2","SUEM840102.polarizability.nf.g3","SUEM840102.second.struct.nf.g1","SUEM840102.second.struct.nf.g2","SUEM840102.second.struct.nf.g3","SUEM840102.solv.acc.nf.g1","SUEM840102.solv.acc.nf.g2","SUEM840102.solv.acc.nf.g3","GEIM800104.hydrophobicity.nf.g1","GEIM800104.hydrophobicity.nf.g2","GEIM800104.hydrophobicity.nf.g3","GEIM800104.vanderwa.nf.g1","GEIM800104.vanderwa.nf.g2","GEIM800104.vanderwa.nf.g3","GEIM800104.charge.nf.g1","GEIM800104.charge.nf.g2","GEIM800104.charge.nf.g3","GEIM800104.polarity.nf.g1","GEIM800104.polarity.nf.g2","GEIM800104.polarity.nf.g3","GEIM800104.polarizability.nf.g1","GEIM800104.polarizability.nf.g2","GEIM800104.polarizability.nf.g3","GEIM800104.second.struct.nf.g1","GEIM800104.second.struct.nf.g2","GEIM800104.second.struct.nf.g3","GEIM800104.solv.acc.nf.g1","GEIM800104.solv.acc.nf.g2","GEIM800104.solv.acc.nf.g3","CHAM830105.hydrophobicity.nf.g1","CHAM830105.hydrophobicity.nf.g2","CHAM830105.hydrophobicity.nf.g3","CHAM830105.vanderwa.nf.g1","CHAM830105.vanderwa.nf.g2","CHAM830105.vanderwa.nf.g3","CHAM830105.charge.nf.g1","CHAM830105.charge.nf.g2","CHAM830105.charge.nf.g3","CHAM830105.polarity.nf.g1","CHAM830105.polarity.nf.g2","CHAM830105.polarity.nf.g3","CHAM830105.polarizability.nf.g1","CHAM830105.polarizability.nf.g2","CHAM830105.polarizability.nf.g3","CHAM830105.second.struct.nf.g1","CHAM830105.second.struct.nf.g2","CHAM830105.second.struct.nf.g3","CHAM830105.solv.acc.nf.g1","CHAM830105.solv.acc.nf.g2","CHAM830105.solv.acc.nf.g3","ARGP820103.hydrophobicity.nf.g1","ARGP820103.hydrophobicity.nf.g2","ARGP820103.hydrophobicity.nf.g3","ARGP820103.vanderwa.nf.g1","ARGP820103.vanderwa.nf.g2","ARGP820103.vanderwa.nf.g3","ARGP820103.charge.nf.g1","ARGP820103.charge.nf.g2","ARGP820103.charge.nf.g3","ARGP820103.polarity.nf.g1","ARGP820103.polarity.nf.g2","ARGP820103.polarity.nf.g3","ARGP820103.polarizability.nf.g1","ARGP820103.polarizability.nf.g2","ARGP820103.polarizability.nf.g3","ARGP820103.second.struct.nf.g1","ARGP820103.second.struct.nf.g2","ARGP820103.second.struct.nf.g3","ARGP820103.solv.acc.nf.g1","ARGP820103.solv.acc.nf.g2","ARGP820103.solv.acc.nf.g3","QIAN880128.hydrophobicity.nf.g1","QIAN880128.hydrophobicity.nf.g2","QIAN880128.hydrophobicity.nf.g3","QIAN880128.vanderwa.nf.g1","QIAN880128.vanderwa.nf.g2","QIAN880128.vanderwa.nf.g3","QIAN880128.charge.nf.g1","QIAN880128.charge.nf.g2","QIAN880128.charge.nf.g3","QIAN880128.polarity.nf.g1","QIAN880128.polarity.nf.g2","QIAN880128.polarity.nf.g3","QIAN880128.polarizability.nf.g1","QIAN880128.polarizability.nf.g2","QIAN880128.polarizability.nf.g3","QIAN880128.second.struct.nf.g1","QIAN880128.second.struct.nf.g2","QIAN880128.second.struct.nf.g3","QIAN880128.solv.acc.nf.g1","QIAN880128.solv.acc.nf.g2","QIAN880128.solv.acc.nf.g3","ARGP820102.hydrophobicity.nf.g1")

output.columns <- c(output.columns,"ARGP820102.hydrophobicity.nf.g2","ARGP820102.hydrophobicity.nf.g3","ARGP820102.vanderwa.nf.g1","ARGP820102.vanderwa.nf.g2","ARGP820102.vanderwa.nf.g3","ARGP820102.charge.nf.g1","ARGP820102.charge.nf.g2","ARGP820102.charge.nf.g3","ARGP820102.polarity.nf.g1","ARGP820102.polarity.nf.g2","ARGP820102.polarity.nf.g3","ARGP820102.polarizability.nf.g1","ARGP820102.polarizability.nf.g2","ARGP820102.polarizability.nf.g3","ARGP820102.second.struct.nf.g1","ARGP820102.second.struct.nf.g2","ARGP820102.second.struct.nf.g3","ARGP820102.solv.acc.nf.g1","ARGP820102.solv.acc.nf.g2","ARGP820102.solv.acc.nf.g3","COSI940101.hydrophobicity.nf.g1","COSI940101.hydrophobicity.nf.g2","COSI940101.hydrophobicity.nf.g3","COSI940101.vanderwa.nf.g1","COSI940101.vanderwa.nf.g2","COSI940101.vanderwa.nf.g3","COSI940101.charge.nf.g1","COSI940101.charge.nf.g2","COSI940101.charge.nf.g3","COSI940101.polarity.nf.g1","COSI940101.polarity.nf.g2","COSI940101.polarity.nf.g3","COSI940101.polarizability.nf.g1","COSI940101.polarizability.nf.g2","COSI940101.polarizability.nf.g3","COSI940101.second.struct.nf.g1","COSI940101.second.struct.nf.g2","COSI940101.second.struct.nf.g3","COSI940101.solv.acc.nf.g1","COSI940101.solv.acc.nf.g2")

output.columns <- c(output.columns,"COSI940101.solv.acc.nf.g3","PALJ810109.hydrophobicity.nf.g1","PALJ810109.hydrophobicity.nf.g2","PALJ810109.hydrophobicity.nf.g3","PALJ810109.vanderwa.nf.g1","PALJ810109.vanderwa.nf.g2","PALJ810109.vanderwa.nf.g3","PALJ810109.charge.nf.g1","PALJ810109.charge.nf.g2","PALJ810109.charge.nf.g3","PALJ810109.polarity.nf.g1","PALJ810109.polarity.nf.g2","PALJ810109.polarity.nf.g3","PALJ810109.polarizability.nf.g1","PALJ810109.polarizability.nf.g2","PALJ810109.polarizability.nf.g3","PALJ810109.second.struct.nf.g1","PALJ810109.second.struct.nf.g2","PALJ810109.second.struct.nf.g3","PALJ810109.solv.acc.nf.g1","PALJ810109.solv.acc.nf.g2","PALJ810109.solv.acc.nf.g3")

# df.input.test <- amp.ind30[1:5,]
impIndexGroup <- c("HUTJ700101","SUEM840102","GEIM800104","CHAM830105","ARGP820103", "QIAN880128","ARGP820102","COSI940101","PALJ810109")
#  df <- amp.ind30[1:5,]
# ind <- impIndexGroup

physico.group <- c('hydrophobicity.g1','hydrophobicity.g2','hydrophobicity.g3','vanderwa.g1','vanderwa.g2','vanderwa.g3','polarity.g1','polarity.g2','polarity.g3','polarizability.g1','polarizability.g2','polarizability.g3','charge.g1','charge.g2','charge.g3','second.struct.g1','second.struct.g2','second.struct.g3','solv.acc.g1','solv.acc.g2','solv.acc.g3')


generate_impFeatures <- function(df, ind){
  
  df.rs <- data.frame(matrix(ncol = 189, nrow = nrow(df)))#change to df
  colnames(df.rs) <- output.columns
  
  for(i in 1:nrow(df)){#O(N)
    seq_length <- df[i,]$seq_length#O(1)
    seq_name <- df[i,]$seq_name#O(1)
    
    hydro.g1.freq <- data.frame(aa = hydrophobicity.g1, qty=str_count(seq_name, hydrophobicity.g1))#O(50)
    hydro.g1.aaindex <- df.aaindex[df.aaindex$aa %in% hydrophobicity.g1, c('aa', ind)]#O(k1)
    hydro.g1.aaindex <- merge(hydro.g1.freq, hydro.g1.aaindex, by = 'aa')#O(k2)
    hydro.g1.aaindex <- hydro.g1.aaindex[hydro.g1.aaindex$qty>0,c('qty', ind)]#O(k3)
    hydro.g1.sum <- sum(hydro.g1.aaindex$qty)#O(k4)
    for(a in 1:9){#O(9)
      df.rs[i,(a-1)*21+1] <- hydro.g1.sum*mean(hydro.g1.aaindex[,ind[a]])/seq_length
    }
    
    hydro.g2.freq <- data.frame(aa = hydrophobicity.g2, qty=str_count(seq_name, hydrophobicity.g2))
    hydro.g2.aaindex <- df.aaindex[df.aaindex$aa %in% hydrophobicity.g2, c('aa', ind)]
    hydro.g2.aaindex <- merge(hydro.g2.freq, hydro.g2.aaindex, by = 'aa')
    hydro.g2.aaindex <- hydro.g2.aaindex[hydro.g2.aaindex$qty>0,c('qty',ind)]
    hydro.g2.sum <- sum(hydro.g2.aaindex$qty)
    for(a in 1:9){
      df.rs[i,(a-1)*21+2] <- hydro.g2.sum*mean(hydro.g2.aaindex[,ind[a]])/seq_length
    }
    
    hydro.g3.freq <- data.frame(aa = hydrophobicity.g3, qty=str_count(seq_name, hydrophobicity.g3))
    hydro.g3.aaindex <- df.aaindex[df.aaindex$aa %in% hydrophobicity.g3, c('aa', ind)]
    hydro.g3.aaindex <- merge(hydro.g3.freq, hydro.g3.aaindex, by = 'aa')
    hydro.g3.aaindex <- hydro.g3.aaindex[hydro.g3.aaindex$qty>0,c('qty',ind)]
    hydro.g3.sum <- sum(hydro.g3.aaindex$qty)
    for(a in 1:9){
      df.rs[i,(a-1)*21+3] <- hydro.g3.sum*mean(hydro.g3.aaindex[,ind[a]])/seq_length
    }
    
    vande.g1.freq <- data.frame(aa = vanderwa.g1, qty=str_count(seq_name, vanderwa.g1))
    vande.g1.aaindex <- df.aaindex[df.aaindex$aa %in% vanderwa.g1, c('aa', ind)]
    vande.g1.aaindex <- merge(vande.g1.freq, vande.g1.aaindex, by = 'aa')
    vande.g1.aaindex <- vande.g1.aaindex[vande.g1.aaindex$qty>0,c('qty',ind)]
    vande.g1.sum <- sum(vande.g1.aaindex$qty)
    for(a in 1:9){
      df.rs[i,(a-1)*21+4] <- vande.g1.sum*mean(vande.g1.aaindex[,ind[a]])/seq_length
    }
    
    vande.g2.freq <- data.frame(aa = vanderwa.g2, qty=str_count(seq_name, vanderwa.g2))
    vande.g2.aaindex <- df.aaindex[df.aaindex$aa %in% vanderwa.g2, c('aa', ind)]
    vande.g2.aaindex <- merge(vande.g2.freq, vande.g2.aaindex, by = 'aa')
    vande.g2.aaindex <- vande.g2.aaindex[vande.g2.aaindex$qty>0,c('qty',ind)]
    vande.g2.sum <- sum(vande.g2.aaindex$qty)
    for(a in 1:9){
      df.rs[i,(a-1)*21+5] <- vande.g2.sum*mean(vande.g2.aaindex[,ind[a]])/seq_length
    }
    
    vande.g3.freq <- data.frame(aa = vanderwa.g3, qty=str_count(seq_name, vanderwa.g3))
    vande.g3.aaindex <- df.aaindex[df.aaindex$aa %in% vanderwa.g3, c('aa', ind)]
    vande.g3.aaindex <- merge(vande.g3.freq, vande.g3.aaindex, by = 'aa')
    vande.g3.aaindex <- vande.g3.aaindex[vande.g3.aaindex$qty>0,c('qty',ind)]
    vande.g3.sum <- sum(vande.g3.aaindex$qty)
    for(a in 1:9){
      df.rs[i,(a-1)*21+6] <- vande.g3.sum*mean(vande.g3.aaindex[,ind[a]])/seq_length
    }
    
    charg.g1.freq <- data.frame(aa = charge.g1, qty=str_count(seq_name, charge.g1))
    charg.g1.aaindex <- df.aaindex[df.aaindex$aa %in% charge.g1, c('aa', ind)]
    charg.g1.aaindex <- merge(charg.g1.freq, charg.g1.aaindex, by = 'aa')
    charg.g1.aaindex <- charg.g1.aaindex[charg.g1.aaindex$qty>0,c('qty',ind)]
    charg.g1.sum <- sum(charg.g1.aaindex$qty)
    for(a in 1:9){
      df.rs[i,(a-1)*21+7] <- charg.g1.sum*mean(charg.g1.aaindex[,ind[a]])/seq_length
    }
    
    charg.g2.freq <- data.frame(aa = charge.g2, qty=str_count(seq_name, charge.g2))
    charg.g2.aaindex <- df.aaindex[df.aaindex$aa %in% charge.g2, c('aa', ind)]
    charg.g2.aaindex <- merge(charg.g2.freq, charg.g2.aaindex, by = 'aa')
    charg.g2.aaindex <- charg.g2.aaindex[charg.g2.aaindex$qty>0,c('qty',ind)]
    charg.g2.sum <- sum(charg.g2.aaindex$qty)
    for(a in 1:9){
      df.rs[i,(a-1)*21+8] <- charg.g2.sum*mean(charg.g2.aaindex[,ind[a]])/seq_length
    }
    
    charg.g3.freq <- data.frame(aa = charge.g3, qty=str_count(seq_name, charge.g3))
    charg.g3.aaindex <- df.aaindex[df.aaindex$aa %in% charge.g3, c('aa', ind)]
    charg.g3.aaindex <- merge(charg.g3.freq, charg.g3.aaindex, by = 'aa')
    charg.g3.aaindex <- charg.g3.aaindex[charg.g3.aaindex$qty>0,c('qty',ind)]
    charg.g3.sum <- sum(charg.g3.aaindex$qty)
    for(a in 1:9){
      df.rs[i,(a-1)*21+9] <- charg.g3.sum*mean(charg.g3.aaindex[,ind[a]])/seq_length
    }
    
    polty.g1.freq <- data.frame(aa = polarity.g1, qty=str_count(seq_name, polarity.g1))
    polty.g1.aaindex <- df.aaindex[df.aaindex$aa %in% polarity.g1, c('aa', ind)]
    polty.g1.aaindex <- merge(polty.g1.freq, polty.g1.aaindex, by = 'aa')
    polty.g1.aaindex <- polty.g1.aaindex[polty.g1.aaindex$qty>0,c('qty',ind)]
    polty.g1.sum <- sum(polty.g1.aaindex$qty)
    for(a in 1:9){
      df.rs[i,(a-1)*21+10] <- polty.g1.sum*mean(polty.g1.aaindex[,ind[a]])/seq_length
    }
    
    polty.g2.freq <- data.frame(aa = polarity.g2, qty=str_count(seq_name, polarity.g2))
    polty.g2.aaindex <- df.aaindex[df.aaindex$aa %in% polarity.g2, c('aa', ind)]
    polty.g2.aaindex <- merge(polty.g2.freq, polty.g2.aaindex, by = 'aa')
    polty.g2.aaindex <- polty.g2.aaindex[polty.g2.aaindex$qty>0,c('qty',ind)]
    polty.g2.sum <- sum(polty.g2.aaindex$qty)
    for(a in 1:9){
      df.rs[i,(a-1)*21+11] <- polty.g2.sum*mean(polty.g2.aaindex[,ind[a]])/seq_length
    }
    
    polty.g3.freq <- data.frame(aa = polarity.g3, qty=str_count(seq_name, polarity.g3))
    polty.g3.aaindex <- df.aaindex[df.aaindex$aa %in% polarity.g3, c('aa', ind)]
    polty.g3.aaindex <- merge(polty.g3.freq, polty.g3.aaindex, by = 'aa')
    polty.g3.aaindex <- polty.g3.aaindex[polty.g3.aaindex$qty>0,c('qty',ind)]
    polty.g3.sum <- sum(polty.g3.aaindex$qty)
    for(a in 1:9){
      df.rs[i,(a-1)*21+12] <- polty.g3.sum*mean(polty.g3.aaindex[,ind[a]])/seq_length
    }
    
    polzt.g1.freq <- data.frame(aa = polarizability.g1, qty=str_count(seq_name, polarizability.g1))
    polzt.g1.aaindex <- df.aaindex[df.aaindex$aa %in% polarizability.g1, c('aa', ind)]
    polzt.g1.aaindex <- merge(polzt.g1.freq, polzt.g1.aaindex, by = 'aa')
    polzt.g1.aaindex <- polzt.g1.aaindex[polzt.g1.aaindex$qty>0,c('qty',ind)]
    polzt.g1.sum <- sum(polzt.g1.aaindex$qty)
    for(a in 1:9){
      df.rs[i,(a-1)*21+13] <- polzt.g1.sum*mean(polzt.g1.aaindex[,ind[a]])/seq_length
    }
    
    polzt.g2.freq <- data.frame(aa = polarizability.g2, qty=str_count(seq_name, polarizability.g2))
    polzt.g2.aaindex <- df.aaindex[df.aaindex$aa %in% polarizability.g2, c('aa', ind)]
    polzt.g2.aaindex <- merge(polzt.g2.freq, polzt.g2.aaindex, by = 'aa')
    polzt.g2.aaindex <- polzt.g2.aaindex[polzt.g2.aaindex$qty>0,c('qty',ind)]
    polzt.g2.sum <- sum(polzt.g2.aaindex$qty)
    for(a in 1:9){
      df.rs[i,(a-1)*21+14] <- polzt.g2.sum*mean(polzt.g2.aaindex[,ind[a]])/seq_length
    }
    
    polzt.g3.freq <- data.frame(aa = polarizability.g3, qty=str_count(seq_name, polarizability.g3))
    polzt.g3.aaindex <- df.aaindex[df.aaindex$aa %in% polarizability.g3, c('aa', ind)]
    polzt.g3.aaindex <- merge(polzt.g3.freq, polzt.g3.aaindex, by = 'aa')
    polzt.g3.aaindex <- polzt.g3.aaindex[polzt.g3.aaindex$qty>0,c('qty',ind)]
    polzt.g3.sum <- sum(polzt.g3.aaindex$qty)
    for(a in 1:9){
      df.rs[i,(a-1)*21+15] <- polzt.g3.sum*mean(polzt.g3.aaindex[,ind[a]])/seq_length
    }
    
    secst.g1.freq <- data.frame(aa = second.struct.g1, qty=str_count(seq_name, second.struct.g1))
    secst.g1.aaindex <- df.aaindex[df.aaindex$aa %in% second.struct.g1, c('aa', ind)]
    secst.g1.aaindex <- merge(secst.g1.freq, secst.g1.aaindex, by = 'aa')
    secst.g1.aaindex <- secst.g1.aaindex[secst.g1.aaindex$qty>0,c('qty',ind)]
    secst.g1.sum <- sum(secst.g1.aaindex$qty)
    for(a in 1:9){
      df.rs[i,(a-1)*21+16] <- secst.g1.sum*mean(secst.g1.aaindex[,ind[a]])/seq_length
    }
    
    secst.g2.freq <- data.frame(aa = second.struct.g2, qty=str_count(seq_name, second.struct.g2))
    secst.g2.aaindex <- df.aaindex[df.aaindex$aa %in% second.struct.g2, c('aa', ind)]
    secst.g2.aaindex <- merge(secst.g2.freq, secst.g2.aaindex, by = 'aa')
    secst.g2.aaindex <- secst.g2.aaindex[secst.g2.aaindex$qty>0,c('qty',ind)]
    secst.g2.sum <- sum(secst.g2.aaindex$qty)
    for(a in 1:9){
      df.rs[i,(a-1)*21+17] <- secst.g2.sum*mean(secst.g2.aaindex[,ind[a]])/seq_length
    }
    
    secst.g3.freq <- data.frame(aa = second.struct.g3, qty=str_count(seq_name, second.struct.g3))
    secst.g3.aaindex <- df.aaindex[df.aaindex$aa %in% second.struct.g3, c('aa', ind)]
    secst.g3.aaindex <- merge(secst.g3.freq, secst.g3.aaindex, by = 'aa')
    secst.g3.aaindex <- secst.g3.aaindex[secst.g3.aaindex$qty>0,c('qty',ind)]
    secst.g3.sum <- sum(secst.g3.aaindex$qty)
    for(a in 1:9){
      df.rs[i,(a-1)*21+18] <- secst.g3.sum*mean(secst.g3.aaindex[,ind[a]])/seq_length
    }
    
    solva.g1.freq <- data.frame(aa = solv.acc.g1, qty=str_count(seq_name, solv.acc.g1))
    solva.g1.aaindex <- df.aaindex[df.aaindex$aa %in% solv.acc.g1, c('aa', ind)]
    solva.g1.aaindex <- merge(solva.g1.freq, solva.g1.aaindex, by = 'aa')
    solva.g1.aaindex <- solva.g1.aaindex[solva.g1.aaindex$qty>0,c('qty',ind)]
    solva.g1.sum <- sum(solva.g1.aaindex$qty)
    for(a in 1:9){
      df.rs[i,(a-1)*21+19] <- solva.g1.sum*mean(solva.g1.aaindex[,ind[a]])/seq_length
    }
    
    solva.g2.freq <- data.frame(aa = solv.acc.g2, qty=str_count(seq_name, solv.acc.g2))
    solva.g2.aaindex <- df.aaindex[df.aaindex$aa %in% solv.acc.g2, c('aa', ind)]
    solva.g2.aaindex <- merge(solva.g2.freq, solva.g2.aaindex, by = 'aa')
    solva.g2.aaindex <- solva.g2.aaindex[solva.g2.aaindex$qty>0,c('qty',ind)]
    solva.g2.sum <- sum(solva.g2.aaindex$qty)
    for(a in 1:9){
      df.rs[i,(a-1)*21+20] <- solva.g2.sum*mean(solva.g2.aaindex[,ind[a]])/seq_length
    }
    
    solva.g3.freq <- data.frame(aa = solv.acc.g3, qty=str_count(seq_name, solv.acc.g3))
    solva.g3.aaindex <- df.aaindex[df.aaindex$aa %in% solv.acc.g3, c('aa', ind)]
    solva.g3.aaindex <- merge(solva.g3.freq, solva.g3.aaindex, by = 'aa')
    solva.g3.aaindex <- solva.g3.aaindex[solva.g3.aaindex$qty>0,c('qty',ind)]
    solva.g3.sum <- sum(solva.g3.aaindex$qty)
    for(a in 1:9){
      df.rs[i,(a-1)*21+21] <- solva.g3.sum*mean(solva.g3.aaindex[,ind[a]])/seq_length
    }
    
  }
  
  
  df.rs[is.na(df.rs)]=0
  #df.rs <- normalizeDF(select(df.rs, -seq_name))
  df.rs
}