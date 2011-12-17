# INFO: Draws chipseq profiles around TSS
# DATE: 28.11.2011.
# AUTH: ABC

# ------------------------------------------- #
# SCRIPT VARIABLES
window.size = 2000

num.of.rand.reg = 5000

peaks.path = "/common/SHARED/vfranke/Fugaku_ChipSeq/Results/Chen_2008_CisgenomePeaks/Peaks/"

outpath="/home/guests/abcsan/SubstituteHistones/Results/DistributionsAroundChenTF"

input.file = "/common/USERS/vfranke/Work/Fugaku_HistoneModification_20052011/Data/ChipSeq/Mapped/RData/Bwa/StandardLengths/FugakuHistones.samp20.w1.uniqTRUE.Normalized.RData"


lib.path = '/home/members/vfranke/Projects/Fugaku_HistoneModification_20052011/Scripts/MasashiLib/SubstituteHistoneChipSeq'
# ------------------------------------------- #
# loads the functions

source(file.path(lib.path,"ProfileFunctions.R"))
source(file.path(lib.path,"Functions.R"))
library(RColorBrewer)

# ------------------------------------------- #
# loads the data
library("BSgenome.Mmusculus.UCSC.mm9")
seqlen = seqlengths(Mmusculus)
Assigner(input.file, "cov.data")

# reading of the peak files and prepearing the data for the analysis
peak.files= list.files(peaks.path, pattern=".cod$", full.names=T, recursive=T)
peaks = lapply(peak.files, read.table, header=T)
names(peaks) = sub("_peak.cod","",basename(peak.files))
l.peaks = list()
for(i in names(peaks)){
    print(i)
    peaks.tmp = peaks[[i]]
    peaks.tmp$start = peaks.tmp$peak_summit - window.size
    peaks.tmp$end   = peaks.tmp$peak_summit + window.size
    peaks.tmp = peaks.tmp[,c("chromosome","start","end", "max_log2FC")]
    names(peaks.tmp)[1] = "chr"
    levs = cut(peaks.tmp$max_log2FC, breaks=fivenum(peaks.tmp$max_log2FC), include.lowest=T)
    levels(levs) = 1:4
    peaks.tmp$set = levs
    l.tmp = list()
    cat("Removing off chromosomes\n")
    for(chr in unique(peaks.tmp$chr)){
        cat(chr,"\r")
        tmp = peaks.tmp[peaks.tmp$chr == chr,]
        l.tmp[[chr]] = tmp[tmp$end < seqlen[chr],]
    }
    peaks[[i]] = do.call(rbind, l.tmp)
}


sample.names = names(cov.data)
b3.ind = grepl("b3", sample.names) 
sample.names = unique(sub("_.+","", sample.names[!b3.ind]))

sample.palette = c(brewer.pal(3, "Set1")[1:2], "darkgray","lightgray")
for(i in 1:length(sample.names)){


    sample = sample.names[i]
    print(sample)
    cov.sample = cov.data[grepl(paste("^",sample,"_", sep=""), names(cov.data)) & !grepl("b3", names(cov.data))]
    names(cov.sample) = sub("h.+","h", names(cov.sample))
    if(length(cov.sample) != 2)
        stop("The number of samples is not 2")

    sample.outpath = file.path(outpath, sample)
        dir.create(sample.outpath, showWarnings=F)

    for(j in 1:length(peaks)){
        
        # ------------------------------------------------------ #
        # converts the chip tf to a GRanges object and generates random regions
        chip.name = names(peaks)[j]
        print(chip.name)
        chip = peaks[[chip.name]]
        chip.ranges = BedToGRanges(chip, values=T)
        o = countOverlaps(chip.ranges, chip.ranges)
        chip.ranges = chip.ranges[o == 1]
        chip.ranges.s = split(chip.ranges, seqnames(chip.ranges))
        chrs = intersect(names(chip.ranges.s), names(cov.sample[[1]]))
        chip.ranges.s = chip.ranges.s[chrs]
        chip.ranges=Reduce(c, chip.ranges.s)
        random.profiles = CreateRandomRegions(seqlen=seqlen[chrs], window.size=2000, ranges=chip.ranges)
        
        
        # ------------------------------------------------------ #
        # gets the coverage for each tf and removes the 10 highest regions
        chip.profiles = lapply(cov.sample, function(x)Coverage2Profiles(x[chrs], chip.ranges.s)) 
        random.profiles = lapply(cov.sample, function(x)Coverage2Profiles(x[chrs], random.profiles))
        # removes overly enriched regions from the chipseq sample - artefacts
        r = lapply(chip.profiles, rowSums)
        r = lapply(r, function(x)head(order(-x), 10))
        chip.profiles.filt = lapply(names(chip.profiles), function(x)chip.profiles[[x]][-r[[x]],])
        names(chip.profiles.filt) = names(chip.profiles)
        
        # removes overly enridched regions from the control datasets
        r.rand = lapply(random.profiles, function(x)head(order(-rowSums(x)), 10))
        rand.profiles.filt = lapply(names(random.profiles), function(x)random.profiles[[x]][-r.rand[[x]],])
        names(rand.profiles.filt) = paste(names(chip.profiles.filt), "rand", sep=".")
        mat.list = c(chip.profiles.filt, rand.profiles.filt)
        
        DrawProfiles(mat.list=mat.list, 
                     fact.list=list(NULL, NULL, NULL, NULL), 
                     indicator.matrix=NULL, 
                     name=paste(sample, chip.name, "w", window.size, "png", sep="."), 
                     outpath=sample.outpath, 
                     palette=sample.palette, 
                     shift=window.size, 
                     split=FALSE)

        # ------------------------------------------------------ #
        # drawing profiles split on the peak height - main quartilles were used
        fact.list = list(paste(chip.name,values(chip.ranges)$set[-r[[1]]], "0h"), paste(chip.name,values(chip.ranges)$set[-r[[2]]], "72h"), NULL, NULL)
        split.cols = c(rep(brewer.pal(3, "Set1")[1:2], each=4), "darkgray","lightgray")
        a = rep(1:4, each=4)
        b = a
        b[seq(2,length(b), 4)] = 5:8
        b[seq(3,length(b), 4)] = 9
        b[seq(4,length(b), 4)] = 10
        ind.mat = cbind(a,b)
        DrawProfiles(mat.list=mat.list, 
                     fact.list=fact.list, 
                     indicator.matrix=ind.mat, 
                     name=paste(sample, chip.name, "w", window.size, "split", "png", sep="."), 
                     outpath=sample.outpath, 
                     palette=split.cols, 
                     shift=window.size, 
                     split=TRUE)


        
    }


}






