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
    granges.tmp = granges.tmp[as.vector(end(granges.tmp)) < seqlen[i] & start(granges.tmp) > 0]
    tss.granges.s[[i]] = granges.tmp
}


#create filename to decide save name
filename<-names(l.data)
chrs = names(seqlen)
for(n in 1:length(filename)){  #repeat
    
    #initialization of matrix
    anno_dataframe=NULL
    sample.name = filename[n]
    cat("Sample name:", sample.name, "\n")
    
    sample.coverage = l.data[[sample.name]]
    # order the samples so that the coverage and the windows have the same order
    cat("Ordering the samples...\n")
    sample.coverage = sample.coverage[order(names(sample.coverage))]
    windows.s = tss.granges.s[order(names(tss.granges.s))]
    strand = unlist(lapply(windows.s, function(x)as.vector(strand(x))))
    if(! all(names(sample.coverage) == names(windows.s)))
        stop("List names do not match")

    cat("Getting the views...\n")
    v = Views(sample.coverage, ranges(windows.s))
    vm = lapply(v, function(x)viewApply(x, as.vector))
    mat = do.call(rbind, lapply(vm, t))

    cat("Reversing the matrices...\n")
    strand.ind = strand == "-"
    mat[strand.ind,] = t(apply(mat[strand.ind,], 1, rev))
    
    tmp.path="/home/guests/abcsan/SubstituteHistones"
    png(file.path(tmp.path, "tryout.png"), width=2000, height=1200)
        plot(-5000:5000, colMeans(mat), col="darkorange", type="l")
    dev.off()

    
    
    
    
    
    