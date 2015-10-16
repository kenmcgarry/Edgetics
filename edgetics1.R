# edgetics1.r
# Using the igraph package and some data from UKCI-2014 paper.
# Edgetics is about using graph theory to model gene mutations as changes
# to network structure by deleting/adding edges to the network. This is turn
# correpsonds to loss of OR gain of protein function.

library(bipartite)
library(Biostrings)
library(protr)
library(rsnps)
library(seqinr)
library(NCBI2R)
library(Peptides)
library(igraph)
library(snplist)
library(biomaRt) # download biomart typing this in console
                  # source("http://bioconductor.org/biocLite.R")
                  # biocLite("bipartite")
data(BLOSUM50)

# THIS SECTION OF CODE PRINTS FIGURE 1 - IN THE PAPER.
el <- read.csv("C:\\R-files\\disease\\DNJuly.csv",header=TRUE)
#el <- el[,1:2] # we just need disease and genes but only first 20

ig <- graph.data.frame(el,directed=FALSE)
V(ig)$type <- V(ig)$name %in% el[,1] # uses matrix mode i.e. %in%
ig

  
#add.vertex.shape("mytriangle"), clip=vertex.shapes("circle")$clip,plot=mytriangle)
add.vertex.shape("triangle", clip=igraph.shape.noclip,plot=mytriangle)

V(ig)$shape <- ifelse(V(ig)$type == TRUE, "circle", "triangle") # select shapes
V(ig)$color <- ifelse(V(ig)$type  == TRUE, "pink", "lightblue")     # select colours

# This is the simple network plot shown in the introduction, save it as jpeg and then open
# in MS paint - resize the plot using the cursor so the borders are close to the nodes,
# otherwise theres a lot of white space that squashes the diagram.

#png("plot1.png", height=6, width=12, units="in", res=200)
#par(mfrow=c(1, 2))


plot.igraph(ig,
vertex.size=25,
vertex.label.color='black',
layout <- layout.fruchterman.reingold(ig),        
            #layout <- layout.reingold.tilford(ig, circular=T),
            layout=layout.fruchterman.reingold(ig, niter=10000, area=30*vcount(ig)^2,axes=F),
            #vertex.label.dist=0.5,
edge.color = "black"
)

#dev.off()

# You cant change the shape of the nodes using tkplot but 
# still could be useful as you can drag the nodes around manually.
tkplot(ig,
       layout = layout.fruchterman.reingold,
       vertex.label = V(ig)$name,
       vertex.label.color= "black",
       ##vertex.size=nodesize,
       vertex.color=V(ig)$color,
       vertex.shape=V(ig)$shape,
       #edge.arrow.size=0,
       edge.curved=FALSE
)

# --------- now do the edgetics --------------

snpInfo<-  cbind(c(17,17,13,13),
                 c(41211653,41213996,32890026,32890572),
                 c("rs8176273","rs8176265","rs9562605","rs1799943"))
colnames(snpInfo)<-c("chr","pos","rsid")
snpInfo <- as.data.frame(snpInfo)

# ----- First obtain the SNP data --------

snp_mart <- useMart("snp", dataset="hsapiens_snp")
snp_attributes = c("refsnp_id", "chr_name", "chrom_start")
snp_locations = getBM(attributes=snp_attributes, filters="snp_filter", 
                      values=snp_ids, mart=snp_mart)

snp_mart <- useMart("snp", dataset="hsapiens_snp")
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
locations <- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol',
                                'chromosome_name', 'start_position', 
                                'end_position'), 
                   filters ='hgnc_symbol', values = "ctns",mart = ensembl)

snps <- getBM(c('refsnp_id','allele','chrom_start','chr_name',
                'distance_to_transcript','ensembl_type','consequence_type_tv'), 
              filters = c('chr_name','chrom_start','chrom_end'), 
              values = list(17,3539762,3564836), mart = snp_mart)

# There are duplications of the rsXXXXX numbers, just keep unique observations.
# REMEMBER - adding or deleting attributes in the call to getBM() will causes the 
# index numbers to change!!!
x<-snps[!duplicated(snps[,1]),]
#x[1:10,] # display first 10 records

# Now select only those identified with HGMD_MUTATION
y<-x[(x[,2])=='HGMD_MUTATION',]
y[1:10,]
# Now write these snips to a CSV file
write.csv(y,file="C:\\R-files\\edgetics\\HGMD_MUTATIONS.csv")

SNPseq<-getBM(c("refsnp_id","snp"),
      filters="snp_filter",
      values=list(c("rs137999539","rs372649305")),mart=snp_mart, checkFilters=TRUE)

utr5 = getSequence(chromosome = 17, start = 3539762, end = 3564836, type = "entrezgene",
                   seqType = "5utr", mart = ensembl)

# getSNP(chromosome = 17, start = 3539762, end = 3564836, mart = snpmart) # some sort of bug??/

#--------------------- using NCBI2R functions 3/7/2014 USE THIS ONE --------
# After using numerous R packages this one seems the best!

#SNPs <- c("rs116000219", "rs116222406", "rs11297981", "rs1119778")
#temp<-NCBI_snp_query(SNPs) # crashes with many snips
#temp<-NCBI_snp_query("rs116222406") # this works...

# GetSNPsInGenes() will return list of snips associated with CTNS (or 1497 in entrez lingo)
## THIS BIT WORKS 26/03/15
SNPList<- GetSNPsInGenes("1497", batchsize=50, MaxRet = 30000, showurl = FALSE,
               quiet = TRUE, smt = FALSE, sme = FALSE)
SNP_info <-GetSNPInfo(SNPList)
# get snips + flanking sequence
mySNPs<-GetSNPFlankSeq(SNPList)   # need to get this right as DNA sequences are amiss

write.csv(x,file="C:\\R-files\\edgetics\\latex_table.csv")
# --------------- Now convert SNP sequences into proteins ----------------------------
# Using seqinr package

SNPseq <- s2c(mySNPs[2,3])   # conversion of a string into a vector "ABC" becomes "A","B","C"
TheSNP <- translate(seq=SNPseq) # Translate nucleic acid sequences into proteins

CTNSseq <- getUniProt("O60931")
CTNSseq <- strsplit(as.character(CTNSseq),"")[[1]] # breakup string into individual aminos 
dotPlot((CTNSseq), TheSNP)

dotPlot(s2c(mySNPs[1,3]),s2c(tim))

x<-length(mySNPs[,3])

# ------ diagnostics ----------
s1 = readFASTA('C://R-files//edgetics//sequence.fasta')[[1]]
s1c <- s2c(s1)
x <- translate(seq=s1c)

# --- seq align from protr package : need to break down into individual bases------
s1<-gsub(", ","",toString(TheSNP))
seqalign = twoSeqSim(CTNSseq, CTNSseq)
summary(seqalign)
print(seqalign@score)


# This is the ORF of CTNS
tim <-"ATGATAAGGAATTGGCTGACTATTTTTATCCTTTTTCCCCTGAAGCTCGTAGAGAAATGTGAGTCAAGCG
TCAGCCTCACTGTTCCTCCTGTCGTAAAGCTGGAGAACGGCAGCTCGACCAACGTCAGCCTCACCCTGCG
GCCACCATTAAATGCAACCCTGGTGATCACTTTTGAAATCACATTTCGTTCCAAAAATATTACTATCCTT
GAGCTCCCCGATGAAGTTGTGGTGCCTCCTGGAGTGACAAACTCCTCTTTTCAAGTGACATCTCAAAATG
TTGGACAACTTACTGTTTATCTACATGGAAATCACTCCAATCAGACCGGCCCGAGGATACGCTTTCTTGT
GATCCGCAGCAGCGCCATTAGCATCATAAACCAGGTGATTGGCTGGATCTACTTTGTGGCCTGGTCCATC
TCCTTCTACCCTCAGGTGATCATGAATTGGAGGCGGAAAAGTGTCATTGGTCTGAGCTTCGACTTCGTGG
CTCTGAACCTGACGGGCTTCGTGGCCTACAGTGTATTCAACATCGGCCTCCTCTGGGTGCCCTACATCAA
GGAGCAGTTTCTCCTCAAATACCCCAACGGAGTGAACCCCGTGAACAGCAACGACGTCTTCTTCAGCCTG
CACGCGGTTGTCCTCACGCTGATCATCATCGTGCAGTGCTGCCTGTATGAGCGCGGTGGCCAGCGCGTGT
CCTGGCCTGCCATCGGCTTCCTGGTGCTCGCGTGGCTCTTCGCATTTGTCACCATGATCGTGGCTGCAGT
GGGAGTGACCACGTGGCTGCAGTTTCTCTTCTGCTTCTCCTACATCAAGCTCGCAGTCACGCTGGTCAAG
TATTTTCCACAGGCCTACATGAACTTTTACTACAAAAGCACTGAGGGCTGGAGCATTGGCAACGTGCTCC
TGGACTTCACCGGGGGCAGCTTCAGCCTCCTGCAGATGTTCCTCCAGTCCTACAACAACGACCAGTGGAC
GCTGATCTTCGGAGACCCAACCAAGTTTGGACTCGGGGTCTTCTCCATCGTCTTCGACGTCGTCTTCTTC
ATCCAGCACTTCTGTTTGTACAGAAAGAGACCGGGGCTTCAGGCAGCGCGCACAGGCTCTGGCAGCCGTC
TCAGGCAGGACTGGGCACCAAGCTTGCAGCCGAAGGCCTTGCCCCAAACTACCAGCGTTTCTGCAAGCAG
CTTGAAGGGCTGA"


#------------------- in-script functions ----------------------
# pretty triangle shape for plotting when using igraph

mytriangle <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size <- 1/150 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  for(i in seq(1,dim(coords)[1])){
    xx <- c(
      coords[i,1],
      coords[i,1]-0.8*vertex.size[i],
      coords[i,1]+0.8*vertex.size[i],
      coords[i,1])
    yy <- c(
      coords[i,2]+0.4*vertex.size[i],
      coords[i,2]-0.8*vertex.size[i],
      coords[i,2]-0.8*vertex.size[i],
      coords[i,2]+0.4*vertex.size[i])
    polygon(x=xx,y=yy,col=vertex.color[i],border="black")
  }
}