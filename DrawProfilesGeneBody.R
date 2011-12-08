# INFO: Draws chipseq whole gene profiles
# DATE: 6.12.2011.
# AUTH: ABC

# ------------------------------------------- #
# loads the functions
source("/home/guests/abcsan/SubstituteHistones/Functions/Functions.R")
source("/home/guests/abcsan/SubstituteHistones/Functions/ProfileFunctions.R")
library(genome.name, character.only=T)
library(RColorBrewer)



# ------------------------------------------- #
# SCRIPT VARIABLES
window.size = 2000

min.gene.length = 4000

num.of.rand.reg = 5000

genome.name = "BSgenome.Mmusculus.UCSC.mm9"

outpath="/home/guests/abcsan/SubstituteHistones/Results/WholeBodyProfiles"

pal = brewer.pal(9, "Set1")

random.col = pal[9]



# ------------------------------------------- #
# loads the data
input.file = "/common/SHARED/vfranke/Fugaku_ChipSeq/Results/NormalizedCoverage/FugakuHistones.samp24.uniqTRUEReducedNormalized.bwa.RData"
Assigner(input.file, "l.data")

input.ind=grepl("input", names(l.data))
cov.data = l.data[!input.ind]
# loads the gene annotation
gene.annotation.path = "/common/SHARED/vfranke/Fugaku_ChipSeq/Results/AnnotatedGenes/Ens.genes.annot.fc2.txt"
gene.annotation<-read.delim(gene.annotation.path, header=T, as.is=T) 
# selects the transcripts with the maximum sum of isoforms
gene.annotation = gene.annotation[order(-rowSums(gene.annotation[,c("X0h.rpkm","X72h.rpkm")])),]
gene.annotation = gene.annotation[!duplicated(gene.annotation$ens.gene.id),]


# ------------------------------------------- #
# removing genes that overlap other genes
genes.ranges = BedToGRanges(gene.annotation, values=T)
start(genes.ranges) = start(genes.ranges) - window.size
end(genes.ranges) = end(genes.ranges) + window.size
count = countOverlaps(genes.ranges, genes.ranges)
genes.selected = genes.ranges[count == 1]
genes.selected = genes.selected[as.vector(seqnames(genes.selected)) == "chr1"]

genes.large = genes.selected[width(genes.selected) > min.gene.length]
tss = ExpandGranges(MakeViewPoint(genes.large, viewpoint="start"), upstream=2000, downstream=20000, strand=TRUE)
gend = ExpandGranges(MakeViewPoint(genes.large, viewpoint="end"), downstream=2000, upstream=20000, strand=TRUE)


### removes the random sequences from the genome lengths
seqlen = seqlengths(Mmusculus)
seqlen = seqlen[!grepl("random", names(seqlen))]
seqlen = seqlen[names(seqlen) != "chrM"]

### removes the windows that fell of the chromosomes
tss.ind = OffChromosomeRegions(tss, seqlen)
gend.ind = OffChromosomeRegions(gend, seqlen)
ind = tss.ind & gend.ind

genes.large = genes.large[ind]
tss = tss[ind]
gend = gend[ind]

# ------------------------------------------- #
# creating random regions on the chromosome
rand.reg = CreateRandomRegions(seqlen, 22000, c(genes.large, tss, gend), n=1, k=1, num.of.rand.reg=1500)

filename<-names(cov.data)
for(n in 1:length(filename)){  #repeat
    
    #initialization of matrix
    anno_dataframe=NULL
    sample.name = filename[n]
    
    # selects the input data for the corresponding time point
    sample.timepoint = str_replace(sample.name,"h.+","h")
    sample.timepoint = str_replace(sample.timepoint,"^.+_","")
    #input.coverage = cov.input[str_detect(names(cov.input), sample.timepoint)][[1]]
    cat("Sample name:", sample.name, "\n")
    
    sample.coverage = cov.data[[sample.name]]
    # order the samples so that the coverage and the windows have the same order
    tss.mat = Coverage2Profiles(sample.coverage["chr1"], split(tss, seqnames(tss))["chr1"])
    gend.mat = Coverage2Profiles(sample.coverage["chr1"], split(gend, seqnames(gend))["chr1"])
    random.mat = Coverage2Profiles(sample.coverage, rand.reg)
    
    random.profile = colMeans(random.mat)
   
#    tss.norm.mat = tss.mat - random.profile
#    gend.norm.mat = t(log2(t(gend.mat+1)/(random.profile+1)))
    
    tss.mat.smooth = t(apply(tss.mat, 1,function(x)ScalerLarge(as.vector(x), 2000)))
    gend.mat.smooth = t(apply(gend.mat, 1, function(x)ScalerLarge(as.vector(x), 2000)))
    rand.mat.smooth = t(apply(random.mat, 1, function(x)ScalerLarge(as.vector(x), 2000)))
#    tss.norm.smooth = t(apply(tss.norm.mat, 1, function(x)ScalerLarge(as.vector(x), 2000)))
#    gend.norm.smooth = t(apply(gend.norm.mat, 1, function(x)ScalerLarge(as.vector(x), 2000)))
    
    tss.norm.smooth = t(log10(t(tss.mat.smooth+1)/(colSums(rand.mat.smooth+1))))
    gene.mat = tss.mat
    #  input.mat  = Coverage2Profiles(input.data, tss.granges.s)
    k = 5
    o = unique(c(head(order(-rowSums(tss.mat)), k), head(order(-rowSums(gend.mat)), k)))
    tss.mat.filt  = tss.mat[-o,]
    gend.mat.filt = gend.mat[-o,]
    tss.smooth.filt = tss.mat.smooth[-o,]
    gend.smooth.filt= gend.mat.smooth[-o,]
    tss.norm.smooth.filt = tss.norm.smooth[-o,]
#    gend.norm.smooth.filt = gend.norm.smooth[-o,]
    
    
    mat.list = list(TSS = tss.mat.filt, GEND=gend.mat.filt)
    tss.to.use = genes.large[-o]

    sample.outpath = file.path(outpath, sample.name)
        dir.create(sample.outpath, showWarnings=F)

    # ----------------------------------------------------- #
    #
    library(marray)
    cols.norm = maPalette(low="darkgray", mid="white", high="darkorange2", k=20)
    cols = colorRampPalette(c("lightgray", "darkblue"), bias=2)(50)
   
    # plots the smoothed profile
    tss.mat.smooth.ord = tss.smooth.filt[order(width(tss.to.use)),]
    png(file.path(sample.outpath, paste(sample.name,"TSSheatmap", "png", sep=".")), width=1000, height=1000)
        image(log10(t(tss.mat.smooth.ord)+1), col=cols, useRaster=T)
    dev.off()
    
    tss.mat.smooth.ord = tss.norm.smooth.filt[order(width(tss.to.use)),]
    tss.mat.smooth.ord = t(log10(apply(tss.mat.smooth.ord,1, function(x)(x-min(x))/(max(x)-min(x)))+1))
    png(file.path(sample.outpath, paste(sample.name,"TSSheatmap","norm", "png", sep=".")), width=1000, height=1000)
        image(t(tss.mat.smooth.ord), col=cols.norm, useRaster=T)
    dev.off()

    gend.mat.smooth.ord =tss.norm.smooth.filt[order(values(tss.to.use)),]
    png(file.path(sample.outpath, paste(sample.name,"GENDheatmap", "png", sep=".")), width=1000, height=1000)
        image(log10(t(gend.mat.smooth.ord)+1), col=cols)
    dev.off()
    
    a = rep(1:2, each=2)
    b=a
    a[seq(2, length(a), 2)] = 3
    indicator.matrix=cbind(b,a)

#    mat.list = list(TSS=tss.mat.filt, Gene.End=gend.mat.filt, Random=random.mat)
    mat.list = list(TSS=tss.mat.filt, Gene.End=gend.mat.filt, Random=random.mat)
    fact.list = list(NULL, NULL, NULL)
    sample.palette = c("darkorange","cornflowerblue","darkgray")

    DrawProfiles(mat.list=mat.list, 
             fact.list=fact.list, 
             indicator.matrix=indicator.matrix, 
             name=paste(sample.name, "png", sep="."), 
             outpath=sample.outpath, 
             palette=sample.palette, 
             shift=c(2000, 20000), 
             split=TRUE)

    
    png(file.path(sample.outpath, paste(sample.name, "png", sep=".")), width=2000, height=800)
#        par(mfrow=c(1,2))
        tss.colmeans = colMeans(tss.mat.filt)
        rand.colmeans = colMeans(random.mat)
        ylim = range(c(range( tss.colmeans), range(rand.colmeans)))
        plot(1:length(tss.colmeans), tss.colmeans, col="darkorange",  type="l", ylim=ylim)
        lines(1:length(rand.colmeans),rand.colmeans, col="darkgray")
    dev.off()

    # ----------------------------------------------------- #
    # bimodal
    bimodal = paste(values(tss.to.use)$H3k4me3_0h, values(tss.to.use)$H3k27me3_0h, sep="_")
    bimodal.f = as.factor(bimodal)
    levels(bimodal.f) = c("None","K27","K4","Bimodal")

    a = rep(1:4, each=2)
    b=a
    a[seq(2, length(a), 2)] = 5
    indicator.matrix=cbind(b,a)

    #    mat.list = list(TSS=tss.mat.filt, Gene.End=gend.mat.filt, Random=random.mat)
    mat.list = list(TSS=tss.mat.filt, Random=random.mat)
    fact.list = list(bimodal.f, NULL)
    sample.palette = c(brewer.pal(4, "Set1"), "darkgray")

    DrawProfiles(mat.list=mat.list, 
                 fact.list=fact.list, 
                 indicator.matrix=indicator.matrix, 
                 name=paste(sample.name,"tss", "bimodal","split","png", sep="."), 
                 outpath=sample.outpath, 
                 palette=sample.palette, 
                 shift=2000, 
                 split=TRUE)


        DrawProfiles(mat.list=mat.list, 
                     fact.list=fact.list, 
                     indicator.matrix=NULL, 
                     name=paste(sample.name,"tss", "bimodal","png", sep="."), 
                     outpath=sample.outpath, 
                     palette=sample.palette, 
                     shift=2000, 
                     split=FALSE)

    
}
