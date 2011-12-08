# INFO: Draws chipseq profiles around TSS
# DATE: 28.11.2011.
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

num.of.rand.reg = 5000

genome.name = "BSgenome.Mmusculus.UCSC.mm9"

outpath="/home/guests/abcsan/SubstituteHistones/Results/ProfileDistributions"

pal = brewer.pal(9, "Set1")

random.col = pal[9]



# ------------------------------------------- #
# loads the data
input.file = "/common/SHARED/vfranke/Fugaku_ChipSeq/Results/NormalizedCoverage/FugakuHistones.samp24.uniqTRUEReducedNormalized.bwa.RData"
Assigner(input.file, "l.data")

input.ind=grepl("input", names(l.data))
cov.input = l.data[input.ind]
cov.data = l.data[!input.ind]
# loads the gene annotation
gene.annotation.path = "/common/SHARED/vfranke/Fugaku_ChipSeq/Results/AnnotatedGenes/Ens.genes.annot.fc2.txt"
gene.annotation<-read.delim(gene.annotation.path, header=T, as.is=T) 
# selects the transcripts with the maximum sum of isoforms
gene.annotation = gene.annotation[order(-rowSums(gene.annotation[,c("X0h.rpkm","X72h.rpkm")])),]
gene.annotation = gene.annotation[!duplicated(gene.annotation$ens.gene.id),]

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
    cat(i,"\r")
    granges.tmp = tss.granges.s[[i]]
    granges.tmp = granges.tmp[as.vector(end(granges.tmp)) < seqlen[i] & start(granges.tmp) > 0]
    tss.granges.s[[i]] = granges.tmp
}
tss.granges = Reduce(c, tss.granges.s)


# creating random regions on the chromosome
wins = MakeTillingWindow(seqlen, (window.size*2+1))
wins.overlap = as.matrix(findOverlaps(wins, genes.expand))
wins.non.overlap = wins[-wins.overlap[,1]]
rand.reg = wins.non.overlap[sample(1:length(wins.non.overlap), num.of.rand.reg)]
rand.reg.s = split(rand.reg, seqnames(rand.reg))



#create filename to decide save name
filename<-names(cov.data)
chrs = names(seqlen)
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
    # first we will do everything without subtraction
    #sample.coverage = sample.coverage - input.coverage
    # order the samples so that the coverage and the windows have the same order
    sample.mat = Coverage2Profiles(sample.coverage, tss.granges.s)
    random.mat = Coverage2Profiles(sample.coverage, rand.reg.s)
    #  input.mat  = Coverage2Profiles(input.data, tss.granges.s)
    r = rowSums(sample.mat)
    o = head(order(-r), 3)
    sample.mat.rem = sample.mat[-o,]
    
    sample.outpath = file.path(outpath, sample.name)
        dir.create(sample.outpath, showWarnings=F)
    
    mat.list = list(sample.mat.rem, Random=random.mat)
    tss.to.use = tss.granges[-o]
    

    # ----------------------------------------------------- #
    # plotting of bimodal data
    bimodal = paste(values(tss.to.use)$H3k4me3_0h, values(tss.to.use)$H3k27me3_0h, sep="_")
    bimodal.f = as.factor(bimodal)
    levels(bimodal.f) = c("None","K27","K4","Bimodal")

    fact.list=list(bimodal.f, NULL)


    
    sample.palette=c(pal[1:4], random.col)
    
    DrawProfiles(mat.list=mat.list, 
                  fact.list=fact.list, 
                  indicator.matrix=NULL, 
                  name=paste(sample.name, "bimodal", "png", sep="."), 
                  outpath=sample.outpath, 
                  palette=sample.palette, 
                  shift=2000, 
                  split=FALSE)
    
    
    a = rep(1:4, each=2)
    b=a
    a[seq(2, length(a), 2)] = 5
    indicator.matrix=cbind(b,a)
     
    DrawProfiles(mat.list=mat.list, 
             fact.list=fact.list, 
             indicator.matrix=indicator.matrix, 
             name=paste(sample.name, "bimodal", "split", "png", sep="."), 
             outpath=sample.outpath, 
             palette=sample.palette, 
             shift=2000, 
             split=TRUE)
    


    # ----------------------------------------------------- #
    # separate.biotypes
    biotypes = values(tss.to.use)$biotype
    biotypes.na = is.na(biotypes)
    biotypes = biotypes[! biotypes.na]
    biotypes.ind = biotypes == "polymorphic_pseudogene" | biotypes == "IG_C_gene" | biotypes == "IG_J_gene" | biotypes == "IG_V_gene"
    biotypes = biotypes[!biotypes.ind]

    sample.mat.rem.biotype = sample.mat.rem[!biotypes.na,]
    sample.mat.rem.biotype = sample.mat.rem.biotype[!biotypes.ind,]
    mat.list.bio = list(sample.mat.rem.biotype, Random=random.mat)


    fact.list=list(biotypes, NULL)


    biotype.palette=c(brewer.pal(9, "Set1"), "black")

    DrawProfiles(mat.list=mat.list.bio, 
                 fact.list=fact.list, 
                 indicator.matrix=NULL, 
                 name=paste(sample.name, "biotype", "png", sep="."), 
                 outpath=sample.outpath, 
                 palette=biotype.palette, 
                 shift=2000, 
                 split=FALSE)

    a = rep(1:length(unique(biotypes)), each=2)
    b=a
    a[1:length(a) %% 2 ==0] = 10
    indicator.matrix.bio=cbind(b,a)

    DrawProfiles(mat.list=mat.list.bio, 
             fact.list=fact.list, 
             indicator.matrix=indicator.matrix.bio[1:8,], 
             name=paste(sample.name, "biotype","split","1-4", "png", sep="."), 
             outpath=sample.outpath, 
             palette=biotype.palette, 
             shift=2000, 
             split=TRUE)

    DrawProfiles(mat.list=mat.list.bio, 
             fact.list=fact.list, 
             indicator.matrix=indicator.matrix.bio[9:18,], 
             name=paste(sample.name, "biotype","split","5-9", "png", sep="."), 
             outpath=sample.outpath, 
             palette=biotype.palette, 
             shift=2000, 
             split=TRUE)



    # -------------------------------------------------------- # 
    # just protein coding genes
#    biotypes = values(tss.to.use)$biotype
#    biotypes.na = is.na(biotypes)
#    biotypes = biotypes[! biotypes.na]
#    biotypes.ind = biotypes == "polymorphic_pseudogene" | biotypes == "IG_C_gene" | biotypes == "IG_J_gene" | biotypes == "IG_V_gene"
#    biotypes = biotypes[!biotypes.ind]
#    
#    bimodal = paste(values(tss.to.use)$H3k4me3_0h, values(tss.to.use)$H3k27me3_0h, sep="_")
#    bimodal.f = as.factor(bimodal)
#    levels(bimodal.f) = c("None","K27","K4","Bimodal")
#    bimodal.f.pc = bimodal.f[!biotypes.na]
#    bimodal.f.pc = bimodal.f.pc[biotypes.pc]
#
#    sample.mat.rem.biotype = sample.mat.rem[!biotypes.na,]
#    sample.mat.rem.biotype.pc = sample.mat.rem.biotype[biotypes.pc,]
#    mat.list.bio.pc = list(sample.mat.rem.biotype.pc, Random=random.mat)
#
#    fact.list.pc=list(bimodal.f.pc, NULL)
#    sample.palette.pc=c(cols, random.col)
#
#    DrawProfiles(mat.list=mat.list.bio.pc, 
#                 fact.list=fact.list.pc, 
#                 indicator.matrix=NULL, 
#                 name=paste(sample.name, "prot.cod","bimodal", "png", sep="."), 
#                 outpath=sample.outpath, 
#                 palette=sample.palette.pc, 
#                 shift=2000, 
#                 split=FALSE)

    # -------------------------------------------------------- # 
    # cpg island--histone modification(0h)

    bimodal = paste(values(tss.to.use)$H3k4me3_0h, values(tss.to.use)$H3k27me3_0h, sep="_")
    cpgs = values(tss.to.use)$cpg.isl.prom
    bimodal.f = as.factor(bimodal)
    levels(bimodal.f) = c("None","K27","K4","Bimodal")
    bimodal.cpg = paste(as.character(bimodal.f), cpgs, sep="_")

    bimodal.cpg.fact.list=list(bimodal.cpg, NULL)
    bimodal.cpg.mat.list = list(sample.mat.rem, Random=random.mat)

    cols = c(brewer.pal(5, "Reds")[4:5], 
             brewer.pal(5, "Purples")[4:5], 
             brewer.pal(5, "Greens")[4:5], 
             brewer.pal(5, "Oranges")[4:5])
    bimodal.cpg.palette=c(cols, random.col)

    DrawProfiles(mat.list=bimodal.cpg.mat.list, 
                 fact.list=bimodal.cpg.fact.list, 
                 indicator.matrix=NULL, 
                 name=paste(sample.name, "bimodal", "cpg", "png", sep="."), 
                 outpath=sample.outpath, 
                 palette=bimodal.cpg.palette, 
                 shift=2000, 
                 split=FALSE)


    a = rep(1:8, each=2)
    b=a
    a[1:length(a) %% 2 ==0] = 9
    indicator.matrix=cbind(b,a)
    DrawProfiles(mat.list=bimodal.cpg.mat.list, 
                 fact.list=bimodal.cpg.fact.list, 
                 indicator.matrix=indicator.matrix[1:8,], 
                 name=paste(sample.name, "bimodal", "cpg", "split", "1-4", "png", sep="."), 
                 outpath=sample.outpath, 
                 palette=bimodal.cpg.palette, 
                 shift=2000, 
                 split=TRUE)

    DrawProfiles(mat.list=bimodal.cpg.mat.list, 
                 fact.list=bimodal.cpg.fact.list, 
                 indicator.matrix=indicator.matrix[9:16,], 
                 name=paste(sample.name, "bimodal", "cpg", "split", "5-8", "png", sep="."), 
                 outpath=sample.outpath, 
                 palette=bimodal.cpg.palette, 
                 shift=2000, 
                 split=TRUE)

    expr= values(tss.to.use)$X0h.rpkm
    expr = log10(expr[-o]+1)
    d= split(data.frame(expr), bimodal.cpg)
    d = lapply(d, unlist)
    png(file=file.path(sample.outpath, "bimodal.cpg.expr.boxplot.png"), width=1000, height=800)
        boxplot(d, bimodal.cpg, col=cols, pch=20, cex=.8, varwidth = T)
    dev.off()

    # -------------------------------------------------------- # 
    #TATA CpG　box--histone modification(0h)
    tata = values(tss.to.use)$tata.box
    tata = ifelse(tata =="Yes", "Tata", "No.Tata")
    
    cpgs = values(tss.to.use)$cpg.isl.prom
    cpgs = ifelse(cpgs=="Yes", "cpg", "No.cpg")
        
    bimodal = paste(values(tss.to.use)$H3k4me3_0h, values(tss.to.use)$H3k27me3_0h, sep="_")
    bimodal.f = as.factor(bimodal)
    levels(bimodal.f) = c("None","K27","K4","Bimodal")
    bimodal.tata.cpg = paste(as.character(bimodal.f), tata, cpgs, sep="_")

    bimodal.tata.cpg.fact.list=list(bimodal.tata.cpg, NULL)
    bimodal.tata.cpg.mat.list = list(sample.mat.rem, Random=random.mat)

    cols=brewer.pal(4, "Set1")
    bimodal.tata.cpg.palette=c(cols, random.col)
    
   
    ind = paste(tata, cpgs, sep="_")[-o]
    
    unique.ind = unique(ind)
    for(i in unique.ind){
        print(i)
        ind.tmp = ind==i
        bimodal.tata.cpg.mat.tmp = bimodal.tata.cpg.mat.list
        bimodal.tata.cpg.mat.tmp[[1]] = bimodal.tata.cpg.mat.tmp[[1]][ind.tmp,]
        bimodal.tata.cpg.fact.tmp = bimodal.tata.cpg.fact.list
        bimodal.tata.cpg.fact.tmp[[1]] = bimodal.tata.cpg.fact.tmp[[1]][ind.tmp]

        DrawProfiles(mat.list=bimodal.tata.cpg.mat.tmp, 
                     fact.list=bimodal.tata.cpg.fact.tmp, 
                     indicator.matrix=NULL, 
                     name=paste(sample.name, "bimodal", i,"png", sep="."), 
                     outpath=sample.outpath, 
                     palette=bimodal.tata.cpg.palette, 
                     shift=2000, 
                     split=FALSE)

    }

    
    # -------------------------------------------------------- # 
    # expression stratification
    expr= values(tss.to.use)$X0h.rpkm
    expr = log10(expr+1)

    expr.class = cut(expr, breaks=c(0,1,4,5,6,max(expr)), include.lowest=T)
    nclass=length(unique(expr.class))

    mat.list.expr = list(sample.mat.rem, Random=random.mat)
    fact.list.expr = list(expr.class, NULL)
    expr.palette = c(brewer.pal(nclass, "Set1"), random.col)

    DrawProfiles(mat.list=mat.list.expr, 
             fact.list=fact.list.expr, 
             indicator.matrix=NULL, 
             name=paste(sample.name, "expr", "png", sep="."), 
             outpath=sample.outpath, 
             palette=expr.palette, 
             shift=2000, 
             split=FALSE)

    a = rep(1:nclass, each=2)
    b= a
    b[1:length(b)%%2 == 0] = nclass+1
    ind.mat = cbind(a,b)
    DrawProfiles(mat.list=mat.list.expr, 
                 fact.list=fact.list.expr, 
                 indicator.matrix=ind.mat, 
                 name=paste(sample.name, "expr", "split","png", sep="."), 
                 outpath=sample.outpath, 
                 palette=expr.palette, 
                 shift=2000, 
                 split=TRUE) 
        
    


    # -------------------------------------------------------- # 
    # expression stratification 0h 72h
    expr.outpath = file.path(sample.outpath, "Expression.0h.72h")
        dir.create(expr.outpath, showWarnings=F)
    expr.0h = log10(values(tss.to.use)$X0h.rpkm+1)
    expr.72h = log10(values(tss.to.use)$X72h.rpkm+1)

    expr.class.0h = paste("0h",cut(expr.0h, breaks=c(0,1,4,5,6,max(expr.0h)), include.lowest=T), sep=" ")
    expr.class.72h = paste("72h", cut(expr.72h, breaks=c(0,1,4,5,6,max(expr.72h)), include.lowest=T), sep=" ")
    unique.ind = unique(expr.class.0h)

    
    a = rep(1:length(unique.ind), each=2)
    b = a
    b[seq(2,length(b),2)] = length(unique.ind)+1
    indicator.matrix = cbind(a,b)

    expr.palette = c(brewer.pal(length(unique.ind), "Set1"), random.col)
    for(i in unique.ind){
        print(i)
        ind = expr.class.0h == i
        expr.mat.list = list(sample.mat.rem[ind,], Random=random.mat)
        expr.fact.list = list(paste(i, expr.class.72h[ind], sep="->"), NULL)
        
        nclass=length(unique(expr.fact.list[[1]]))
        a = rep(1:nclass, each=2)
        b = a
        b[seq(2,length(b),2)] = nclass+1
        indicator.matrix = cbind(a,b)

        expr.palette = c(brewer.pal(nclass, "Set1"), random.col)
        DrawProfiles(mat.list=expr.mat.list, 
                    fact.list=expr.fact.list, 
                    indicator.matrix=NULL, 
                    name=paste(sample.name, "expr", i,"png", sep="."), 
                    outpath=expr.outpath, 
                    palette=expr.palette, 
                    shift=2000, 
                    split=FALSE)
        
        DrawProfiles(mat.list=expr.mat.list, 
                     fact.list=expr.fact.list, 
                     indicator.matrix=indicator.matrix, 
                     name=paste(sample.name, "expr", "split", i,"png", sep="."), 
                     outpath=expr.outpath, 
                     palette=expr.palette, 
                     shift=2000, 
                     split=TRUE)

    
    # -------------------------------------------------------- # 
    # bimodal 0h - bimodal 72h
    bimod.outpath = file.path(sample.outpath, "Bimodality.0h.72h")
        dir.create(bimod.outpath, showWarnings=F)
    bimodal.0h = as.factor(paste(values(tss.to.use)$H3k4me3_0h, values(tss.to.use)$H3k27me3_0h, sep="_"))
    levels(bimodal.0h) = paste("0h",c("None","K27","K4","Bimodal"))

    bimodal.72h = as.factor(paste(values(tss.to.use)$H3k4me3_72h, values(tss.to.use)$H3k27me3_72h, sep="_"))
    levels(bimodal.72h) = paste("72h",c("None","K27","K4","Bimodal"))

    unique.ind = unique(bimodal.0h)

    for(i in unique.ind){
        print(i)
        ind = bimodal.0h == i
        bimod.mat.list = list(sample.mat.rem[ind,], Random=random.mat)
        bimod.fact.list = list(paste(i, bimodal.72h[ind], sep="->"), NULL)
            
        nclass=length(unique(bimod.fact.list[[1]]))
        a = rep(1:nclass, each=2)
        b = a
        b[seq(2,length(b),2)] = nclass+1
        indicator.matrix = cbind(a,b)
            
        bimod.palette = c(brewer.pal(nclass, "Set1"), random.col)
        DrawProfiles(mat.list=bimod.mat.list, 
                    fact.list=bimod.fact.list, 
                    indicator.matrix=NULL, 
                    name=paste(sample.name, "bimod", i,"png", sep="."), 
                    outpath=bimod.outpath, 
                    palette=bimod.palette, 
                    shift=2000, 
                    split=FALSE)
        
        DrawProfiles(mat.list=bimod.mat.list, 
                    fact.list=bimod.fact.list, 
                    indicator.matrix=indicator.matrix, 
                    name=paste(sample.name, "bimod", "split", i,"png", sep="."), 
                    outpath=bimod.outpath, 
                    palette=bimod.palette, 
                    shift=2000, 
                    split=TRUE)


    }




cat("Everything went ok! ... bye! \n")
    
    
    