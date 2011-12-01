# INFO: Draws chipseq profiles around TSS
# DATE: 28.11.2011.
# AUTH: ABC

# ------------------------------------------- #
# SCRIPT VARIABLES
window.size = 2000

num.of.rand.reg = 5000

genome.name = "BSgenome.Mmusculus.UCSC.mm9"

outpath="/home/guests/abcsan/SubstituteHistones/Results/ProfileDistributions"


# ------------------------------------------- #
# loads the functions
source("/home/guests/abcsan/SubstituteHistones/Functions/Functions.R")
source("/home/guests/abcsan/SubstituteHistones/Functions/ProfileFunctions.R")
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

### removes granges that overlapp another feature on the genome
tss.granges = BedToGRanges(tss.window, values=T)
genes.expand = RegionExtender(gene.annotation, downstream=window.size, upstream=window.size)
genes.expand = BedToGRanges(genes.expand)
o = as.matrix(findOverlaps(tss.granges, genes.expand, ignore.strand = T))
o.ind = tapply(o[,2], o[,1], length)
tss.granges.non = tss.granges[as.numeric(names(o.ind)[o.ind == 1])]

tss.granges.s = split(tss.granges.non, seqnames(tss.granges.non))

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


# creating random regions on the chromosome
wins = MakeTillingWindow(seqlen, window.size*2)
wins.overlap = as.matrix(findOverlaps(wins, genes.expand))
wins.non.overlap = wins[-wins.overlap[,1]]
rand.reg = wins.non.overlap[sample(1:length(wins.non.overlap), num.of.rand.reg)]
rand.reg.s = split(rand.reg, seqnames(rand.reg))


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
    sample.mat = Coverage2Profiles(sample.coverage, tss.granges.s)
    random.mat = Coverage2Profiles(sample.coverage, rand.reg.s)
    
    
    # plots the cumulative profile
    cat("Drawing profiles...\n")
    png(file.path(outpath, paste(sample.name, "w", window.size, "png", sep=".")), width=2000, height=1200)
        plot(-window.size:window.size, colMeans(mat), col="darkorange", type="l", lwd=2, ylab="mean.coverage", xlab="position", main=sample.name)
    dev.off()
    
    # splits the TSS by valency
    bimodal = unlist(lapply(windows.s, function(x)paste(values(x)$H3k4me3_0h, values(x)$H3k27me3_0h, sep="_")))

    mat.s = split(as.data.frame(mat), bimodal)
    mat.s = lapply(mat.s, colMeans)

    # plots separately profiles for each modality
    cols = c("darkorange","cornflowerblue","firebrick","chartreuse")
    png(file.path(outpath, paste(sample.name, "bivalent", "w", window.size, "png", sep=".")), width=2000, height=2000)
        par(mfrow=c(2,2), cex.main=2.5, cex.lab=2.5, cex.axis=2.5)
        for(i in 1:length(mat.s)){
            print(i)
            plot(-window.size:window.size, mat.s[[i]], col=cols[i], main=names(mat.s)[i], type="l", lwd=2, ylab="mean.coverage", xlab="position")
        }

    dev.off()

    # plots modality profiles on one plot
    cols = c("darkorange","cornflowerblue","firebrick","chartreuse")
    png(file.path(outpath, paste(sample.name, "bivalent", "oneplot", "w", window.size, "png", sep=".")), width=1200, height=1200)
        par(cex.main=2, cex.lab=2, cex.axis=2)
        for(i in 1:length(mat.s)){
            print(i)
                                                 
            if(i == 1){
                                                
                plot(-window.size:window.size, mat.s[[i]], col=cols[i], type="l", lwd=2, ylim=c(0, max(unlist(mat.s))), main=sample.name, ylab="mean.coverage", xlab="position")
                abline(v=0, lwd=2)
            }else{
                lines(-window.size:window.size, mat.s[[i]], col=cols[i], lwd=2, )                                 
            }
        }
        legend("topright", fill=cols, legend=names(mat.s))                                         

    dev.off()

    
}

cat("Everything went ok! ... bye! \n")
    
    
    