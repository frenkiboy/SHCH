# INFO: Draws chipseq profiles around TSS
# DATE: 28.11.2011.
# AUTH: ABC

# ------------------------------------------- #
# SCRIPT VARIABLES
window.size = 5000

genome.name = "BSgenome.Mmusculus.UCSC.mm9"


# ------------------------------------------- #
# loads the functions
source("/home/guests/abcsan/SubstituteHistones/Functions/Functions.R")
library(genome.name, character.only=T)


# ------------------------------------------- #
# loads the data
input.file = "/common/SHARED/vfranke/Fugaku_ChipSeq/Results/NormalizedCoverage/FugakuHistones.samp24.uniqTRUEReducedNormalized.bwa.RData"
Assigner(input.file, "l.data")

# loads the gene annotation
gene.annotation.path = "/common/SHARED/vfranke/Fugaku_ChipSeq/Results/AnnotatedGenes/Ens.genes.annot.fc2.txt"
gene.annotation<-read.delim(gene.annotation.path, header=T, as.is=T) 


# ------------------------------------------- #
# check the sample names 
sub("h.+","h", names(l.data))

tss = StartStrandMaker(gene.annotation)
tss.window = RegionExtender(tss, downstream=window.size, upstream=window.size)
tss.granges = BedToGRanges(tss.window, values=T)
tss.granges.s = split(tss.granges, seqnames(tss.granges))

### removes the random sequences from the genome lengths
seqlen = seqlengths(Mmusculus)
seqlen = seqlen[!grepl("random", names(seqlen))]
seqlen = seqlen[names(seqlen) != "chrM"]

### removes the windows that fell of the chromosomes
for(i in names(seqlen)){
    cat(i,"\n")
    granges.tmp = tss.granges.s[[i]]
    granges.tmp = granges.tmp[as.vector(end(granges.tmp)) < seqlen[i] & start(granges) > 0]
    tss.granges.s[[i]] = granges.tmp
}


#create filename to decide save name
filename<-names(l.data)
chrs = names(seqlen)
for(n in 1:length(filename)){  #repeat
    
    #initialization of matrix
    anno_dataframe=NULL
    sample.coverage = l.data[[filename[n]]]
    # order the samples so that the coverage and the windows have the same order
    sample.coverage = sample.coverage[order(chrs)]
    windows.s = tss.granges.s[order(chrs)]
    
    v = Views(sample.coverage, ranges(windows.s))
    v = lapply(v, function(x)viewApply(x, as.vector))
    mat = do.call(rbind, lapply(v, t))
    
    
    
    
    
    
    
    for(k in 1:21){     #repeat 
        chromosome <-names(l.norm[[n]][k])     #choose chromosome
        chr_anno<-annotation[annotation[ ,1]==chromosome, ]   #select the gene whose chromosome correspond to chromosome of repeat
        for(i in 1:nrow(chr_anno)){
            chr_cov<-Views(l.norm[[n]][[k]], start=chr_anno[i, 2]-5000, end=chr_anno[i, 2]+5000, name=chr_anno[i, 6])  #these names are derived from ensemble transcript ids. they are not duplicated.
            chr_vector<-as.vector(chr_cov[[1]])
            chr_matrix<-as.matrix(chr_vector)
            colnames(chr_matrix)<-chr_anno[i, 6]
            anno_dataframe<-cbind(anno_dataframe, chr_matrix)
    }
    }
        write.table(anno_dataframe, file.path("/home/ymasashi/ABC/R_and_perl_for_ChIPseq_analysis/", paste(filename[n], "TSS+-5kb_annotation.polyA-_tpm.txt", sep=".")), col.names=T, row.names=F, quote=F, sep="\t")   #save TSS+- 5kb coverage as dataframe. automatically decide file name by file.path function.
    gc()
}