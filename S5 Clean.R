#Supplementary Figure 5
#Download sequences as detailed in the readme.  I then loaded them into Geneious and mapped to the reference sequence.
#I exported the alignments as Fasta files to my working directory.
#The following script removes all reads with >10 errors compared to the reference.  It then analyses the remaining reads can calculates whether the mutations are mostly in the first 10 or last 10 of a read which implies that it is sequencing error or DIs.
#Remaining mutations were analysed in Geneious by eye to confirm the mutations.


#Sequencing analysis
library(QuasR)
library(ShortRead)
library(seqinr)
filey=list.files(pattern="fasta") #Or use an appropriate pattern to find your files
for(fn in 1:length(filey))
{
  spl=0
  filez=filey[fn]
  cat(filez,"\n")
  cat(filez,"\n",file="output.txt",append=TRUE)
  #alnss=readDNAStringSet(filez)
  alnrl=readLines(filez,2)[2]
  aln=readDNAMultipleAlignment(filez)
  if(dim(aln)[1]>100000) #I had to put this bit in to stop the program crashing with giant alignments.  It splits it into bits of 100,000.
  {
    alnss=readDNAStringSet(filez)
    spl=1
  }
  if(spl==0)
  {
    gap_positions <- which(strsplit(alnrl, "")[[1]] == "-") #Identify gaps in alignment
    cat("Making the matrix")
    colmask(aln)=IRanges(gap_positions,gap_positions,1)
    mataln=as.matrix(aln)
    rm(aln)
    
    #Count mismatches and remove seqs >10 mismatches
    ref <- mataln[1, ]
    mismatch <- t(apply(mataln, 1, function(row) {
      m <- as.integer(row != ref)
    }))
    #Mismatches are 1. This replaces - with 0.
    binary_mat <- ifelse(mataln == "-", 0, 1)
    mismatch=mismatch*binary_mat
    badseqs=which(rowSums2(mismatch)>10)
    mataln=mataln[-badseqs,]
    binary_mat=binary_mat[-badseqs,]
    mismatch=mismatch[-badseqs,]
    #Find more bad seqs with stuff too close to the end
    last_one <- max.col(binary_mat, ties.method = "last")
    sum_last10 <- mapply(function(i, last) {
      
      start <- max(1, last - 9)             # start of the 10-column window
      sum(mismatch[i, start:last])               # sum within that window
    }, i = seq_len(nrow(mismatch)), last = last_one)
    
    first_one=max.col(binary_mat, ties.method = "first")
    #Find more bad seqs with stuff too close to the start
    sum_first10 <- mapply(function(i, first) {
      actlast1=min(first+9,ncol(mismatch))
      sum(mismatch[i, first:actlast1])               # sum within that window
    }, i = seq_len(nrow(mismatch)),first=first_one)
    badseqs2=union(which(sum_last10>2),which(sum_first10>2))
    mataln=mataln[-badseqs2,]
    binary_mat=binary_mat[-badseqs2,]
    mismatch=mismatch[-badseqs2,]
    last_one=last_one[-badseqs2]
    first_one=first_one[-badseqs2]
    aln_new <- DNAStringSet(apply(mataln, 1, paste0, collapse = ""))
    
    #Make a consensus matrix
    mat=consensusMatrix(aln_new)[1:4,]
    findAmbiguousColumns <- function(mat, threshold = 0.04) {
      apply(mat, 2, function(col) {
        freqs <- col / sum(col)
        sum(freqs > threshold) >= 2
      }) |> which()
    }
    ambig_columns=findAmbiguousColumns(mat,0.04)
    if(length(ambig_columns)>0.5)
    {
      cat("Found some possible columns",ambig_columns,"\n")
      cat("Found some possible columns",ambig_columns,"\n",file="output.txt",append=TRUE)
      #Check those columns
      for(a in 1:length(ambig_columns))
      {
        second_base=order(mat[,ambig_columns[a]])[3]
        if(mataln[1,ambig_columns[a]]==c("A","C","G","T")[second_base])
        {cat("Are we returning to the reference - check this for",ambig_columns[a])
          cat("Are we returning to the reference - check this for",ambig_columns[a],file="output.txt",append=TRUE)
        }
        baddies=which(mismatch[,ambig_columns[a]]==1)
        cat("This percentage of",length(baddies),"sequences in ",ambig_columns[a], "is in the first or last bit-",length(which(first_one[baddies]%in%c((ambig_columns[a]-9):ambig_columns[a])|last_one[baddies]%in%c((ambig_columns[a]+9):ambig_columns[a])))/length(baddies),"\n")  
        cat("This percentage of",length(baddies),"sequences in ",ambig_columns[a], "is in the first or last bit-",length(which(first_one[baddies]%in%c((ambig_columns[a]-9):ambig_columns[a])|last_one[baddies]%in%c((ambig_columns[a]+9):ambig_columns[a])))/length(baddies),"\n",file="output.txt",append=TRUE)  
      }
    }
    if(length(ambig_columns)<0.5)
    {
      cat("Found no ambiguous columns","\n")
      cat("Found no ambiguous columns","\n",file="output.txt",append=TRUE)
    }
    file_name=paste(c(strsplit(filez," ")[[1]][1],"_adj.fasta"),collapse="")
    writeXStringSet(aln_new, file_name)
  }
  if(spl==1)
  {
    bignum=dim(aln)[1]
    splitarray=array(dim=c(ceiling(bignum/100000),2))
    for(aa in 1:ceiling(bignum/100000))
    {
      splitarray[aa,1]=aa*100000-99999
      splitarray[aa,2]=aa*100000
    }
    splitarray[ceiling(bignum/100000),2]=bignum
    for(aa in 1:ceiling(bignum/100000))
    {
      cat("Split, split, split",aa,"\n")
      cat("Split, split, split",aa,"\n",file="output.txt",append=TRUE)
      aln=DNAMultipleAlignment(alnss[c(1,splitarray[aa,1]:splitarray[aa,2]),])
      gap_positions <- which(strsplit(alnrl, "")[[1]] == "-") #Identify gaps in alignment
      cat("Making the matrix")
      colmask(aln)=IRanges(gap_positions,gap_positions,1)
      mataln=as.matrix(aln)
      rm(aln)
      
      #Count mismatches and remove seqs >10 mismatches
      ref <- mataln[1, ]
      mismatch <- t(apply(mataln, 1, function(row) {
        m <- as.integer(row != ref)
      }))
      #Mismatches are 1. This replaces - with 0.
      binary_mat <- ifelse(mataln == "-", 0, 1)
      mismatch=mismatch*binary_mat
      badseqs=which(rowSums2(mismatch)>10)
      mataln=mataln[-badseqs,]
      binary_mat=binary_mat[-badseqs,]
      mismatch=mismatch[-badseqs,]
      #Find more bad seqs with stuff too close to the end
      last_one <- max.col(binary_mat, ties.method = "last")
      sum_last10 <- mapply(function(i, last) {
        
        start <- max(1, last - 9)             # start of the 10-column window
        sum(mismatch[i, start:last])               # sum within that window
      }, i = seq_len(nrow(mismatch)), last = last_one)
      
      first_one=max.col(binary_mat, ties.method = "first")
      #Find more bad seqs with stuff too close to the start
      sum_first10 <- mapply(function(i, first) {
        actlast1=min(first+9,ncol(mismatch))
        sum(mismatch[i, first:actlast1])               # sum within that window
      }, i = seq_len(nrow(mismatch)),first=first_one)
      badseqs2=union(which(sum_last10>2),which(sum_first10>2))
      mataln=mataln[-badseqs2,]
      binary_mat=binary_mat[-badseqs2,]
      mismatch=mismatch[-badseqs2,]
      last_one=last_one[-badseqs2]
      first_one=first_one[-badseqs2]
      aln_new <- DNAStringSet(apply(mataln, 1, paste0, collapse = ""))
      
      #Make a consensus matrix
      mat=consensusMatrix(aln_new)[1:4,]
      findAmbiguousColumns <- function(mat, threshold = 0.04) {
        apply(mat, 2, function(col) {
          freqs <- col / sum(col)
          sum(freqs > threshold) >= 2
        }) |> which()
      }
      ambig_columns=findAmbiguousColumns(mat,0.04)
      if(length(ambig_columns)>0.5)
      {
        cat("Found some possible columns",ambig_columns,"\n")
        cat("Found some possible columns",ambig_columns,"\n",file="output.txt",append=TRUE)
        #Check those columns
        for(a in 1:length(ambig_columns))
        {
          second_base=order(mat[,ambig_columns[a]])[3]
          if(mataln[1,ambig_columns[a]]==c("A","C","G","T")[second_base])
          {cat("Are we returning to the reference - check this for",ambig_columns[a])
            cat("Are we returning to the reference - check this for",ambig_columns[a],file="output.txt",append=TRUE)
          }
          baddies=which(mismatch[,ambig_columns[a]]==1)
          cat("This percentage of",length(baddies),"sequences in ",ambig_columns[a], "is in the first or last bit-",length(which(first_one[baddies]%in%c((ambig_columns[a]-9):ambig_columns[a])|last_one[baddies]%in%c((ambig_columns[a]+9):ambig_columns[a])))/length(baddies),"\n")  
          cat("This percentage of",length(baddies),"sequences in ",ambig_columns[a], "is in the first or last bit-",length(which(first_one[baddies]%in%c((ambig_columns[a]-9):ambig_columns[a])|last_one[baddies]%in%c((ambig_columns[a]+9):ambig_columns[a])))/length(baddies),"\n",file="output.txt",append=TRUE)  
        }
      }
      if(length(ambig_columns)<0.5)
      {
        cat("Found no ambiguous columns","\n")
        cat("Found no ambiguous columns","\n",file="output.txt",append=TRUE)
      }
      file_name=paste(c(strsplit(filez," ")[[1]][1],aa,"_adj.fasta"),collapse="")
      writeXStringSet(aln_new, file_name)
      
    }
  }
}
#Make a nice graph (Supplementary Figure S5)
#In truth, this makes 3 graphs which I cut and pasted to make the actual figure because I couldn't be bothered to work out how to work it out using ggplot.
samp <- c(
  "Cow/Texas STOCK", 
  "Cow/Texas hNAEC #164 72h (A)", "Cow/Texas hNAEC #164 72h (B)", "Cow/Texas hNAEC #164 72h (C)", 
  "Cow/Texas hNAEC #226 72h (A)", "Cow/Texas hNAEC #226 72h (B)", "Cow/Texas hNAEC #226 72h (C)", 
  "Cow/Texas PB2-D740N STOCK", 
  "Cow/Texas PB2-D740N hNAEC #164 72h (A)", "Cow/Texas PB2-D740N hNAEC #164 72h (B)", "Cow/Texas PB2-D740N hNAEC #164 72h (C)", 
  "Cow/Texas PB2-D740N hNAEC #226 72h (A)", "Cow/Texas PB2-D740N hNAEC #226 72h (B)", "Cow/Texas PB2-D740N hNAEC #226 72h (C)", 
  "Cow/Texas avianised polymerase STOCK", 
  "Cow/Texas avianised polymerase hNAEC #164 72h (A)", "Cow/Texas avianised polymerase hNAEC #164 72h (B)", "Cow/Texas avianised polymerase hNAEC #164 72h (C)", 
  "Cow/Texas avianised polymerase hNAEC #226 72h (A)", "Cow/Texas avianised polymerase hNAEC #226 72h (B)", "Cow/Texas avianised polymerase hNAEC #226 72h (C)"
)

# Define the shape sequence (one per population, in same order)
shape_seq <- c(19,15,17,18,15,17,18,19,15,17,18,15,17,18,19,15,17,18,15,17,18)
sample_colours <- scales::hue_pal()(21)[c(1,3,3,3,6,6,6,8,10,10,10,13,13,13,15,17,17,17,20,20,20)]
names(sample_colours) <- samp
HAdata=data.frame(
  x =c(32	,219,	235,235,235	,384,	385,	425,	553),y=c(11,5,15,31,37,20,8,5,6),aalabel=c("V32I","S219L","","","V235I","E384K","S385Y","A448E","L553S"),sample = factor(c(samp[c(2,11,11,13,14,21,2,19,13)]), levels = samp))
PAdata=data.frame(
  x =c(206,	260	,290,	418,	431,	448	,636),y=c(8,7.5,12,10,4,22,6),aalabel=c("E206G","L260F","L290V","T418A","D431N","A448E","V636A"),sample = factor(c(samp[c(8,7,5,10,20,21,14)]), 
                                                                                                                                                     levels = samp))
PB2data <- data.frame(
  x = c(102, 478, 652),
  y = c(4, 9, 4),aalabel=c("N102D","I478M","N652D"),
  sample = factor(c(samp[c(5,21,16)]), 
                  levels = samp)  # include all 21 in factor levels
)
PB2data$sample <- factor(PB2data$sample, levels = samp)
# Create a named shape vector for scale_shape_manual
shape_map <- setNames(shape_seq, samp)
levels_samples <- levels(PAdata$sample)

# create a named vector mapping each sample level to the corresponding shape
shape_map <- setNames(shape_seq, levels_samples)




HAdata$sample <- factor(HAdata$sample, levels = samp)

# --- Plot ---
ggplot(HAdata, aes(x = x, y = y, colour = sample, shape = sample, label = aalabel)) +
  theme_minimal() +
  geom_point(size = 3) +
  geom_segment(aes(xend = x, y = 0, yend = y)) +
  geom_text(vjust = -1.5, size = 3) +
  geom_hline(yintercept = 4, linetype = "dashed", colour = "black") +
  scale_x_continuous(limits = c(1, 568)) +
  scale_y_continuous(limits = c(0, 100)) +
  scale_colour_manual(values = sample_colours, drop = FALSE) +
  scale_shape_manual(values = shape_map, drop = FALSE) +
  labs(
    x = "HA amino acid position (immature H5 numbering)",
    y = "% Frequency"
  )

ggplot(PAdata, aes(x = x, y = y, colour = sample, shape = sample, label = aalabel)) +
  theme_minimal() +
  geom_point(size = 3) +
  geom_segment(aes(xend = x, y = 0, yend = y)) +
  geom_text(vjust = -1.5, size = 3) +
  geom_hline(yintercept = 4, linetype = "dashed", colour = "black") +
  scale_x_continuous(limits = c(1, 717)) +
  scale_y_continuous(limits = c(0, 100)) +
  scale_colour_manual(values = sample_colours, drop = FALSE) +
  scale_shape_manual(values = shape_map, drop = FALSE) +
  labs(
    x = "PA amino acid position",
    y = "% Frequency"
  )

ggplot(PB2data, aes(x = x, y = y, colour = sample, shape = sample, label = aalabel)) +
  theme_minimal() +
  geom_point(size = 3) +
  geom_segment(aes(xend = x, y = 0, yend = y)) +
  geom_text(vjust = -1.5, size = 3) +
  geom_hline(yintercept = 4, linetype = "dashed", colour = "black") +
  scale_x_continuous(limits = c(1, 760)) +
  scale_y_continuous(limits = c(0, 100)) +
  scale_colour_manual(values = sample_colours, drop = FALSE) +
  scale_shape_manual(values = shape_map, drop = FALSE) +
  labs(
    x = "PB2 amino acid position",
    y = "% Frequency"
  )

