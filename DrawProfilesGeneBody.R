#!/common/software/shared/R-2.14.0/bin/R
# INFO: Draws chipseq whole gene profiles
# DATE: 6.12.2011.
# AUTH: ABC

# ------------------------------------------- #
# loads the functions
functions.path = '/home/members/vfranke/Projects/Fugaku_HistoneModification_20052011/Scripts/MasashiLib/SubstituteHistoneChipSeq'
source(file.path(functions.path, "Functions.R"))
source(file.path(functions.path, "ProfileFunctions.R"))
library(RColorBrewer)



# ------------------------------------------- #
# SCRIPT VARIABLES
window.size = 2000

min.gene.length = 4000

num.of.rand.reg = 5000

genome.name = "BSgenome.Mmusculus.UCSC.mm9"

outpath="/common/SHARED/vfranke/Fugaku_ChipSeq/Results/GeneBodyProfiles"

pal = brewer.pal(9, "Set1")

random.col = pal[9]

input.file = "/common/SHARED/vfranke/Fugaku_ChipSeq/Results/NormalizedCoverage/FugakuHistones.samp24.uniqTRUEReducedNormalized.bwa.RData"

gene.annotation.path = "/common/SHARED/vfranke/Fugaku_ChipSeq/Results/AnnotatedGenes/Ens.genes.annot.fc2.txt"

# ------------------------------------------- #
# loads the data
Assigner(input.file, "cov.data")
genome = GenomeLoader(genome.name)

# loads the gene annotation
gene.annotation = read.delim(gene.annotation.path, header=T, as.is=T) 
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
genes.selected = genes.selected[as.vector(seqnames(genes.selected)) %in% c("chr1")]

genes.large = genes.selected[width(genes.selected) > min.gene.length]
tss = ExpandGranges(MakeViewPoint(genes.large, viewpoint="start"), upstream=2000, downstream=20000, strand=TRUE)
gend = ExpandGranges(MakeViewPoint(genes.large, viewpoint="end"), downstream=2000, upstream=20000, strand=TRUE)


### removes the random sequences from the genome lengths
seqlen = seqlengths(genome)
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

filename = names(cov.data)
for(n in 1:length(filename)){
    
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
    
	tss.mat.smooth = t(apply(tss.mat, 1,function(x)ScalerLarge(as.vector(x), 2000)))
    gend.mat.smooth = t(apply(gend.mat, 1, function(x)ScalerLarge(as.vector(x), 2000)))
    rand.mat.smooth = t(apply(random.mat, 1, function(x)ScalerLarge(as.vector(x), 2000)))
    
	tss.norm.smooth = t(log10(t(tss.mat.smooth+1)/(colMeans(rand.mat.smooth+1))))
    gend.norm.smooth = t(log10(t(gend.mat.smooth+1)/(colMeans(rand.mat.smooth+1)))) 
    # removes very high coverage regions
    k = 5
    o = unique(c(head(order(-rowSums(tss.mat)), k), head(order(-rowSums(gend.mat)), k)))
    tss.mat.filt  = tss.mat[-o,]
    gend.mat.filt = gend.mat[-o,]
    tss.smooth.filt = tss.mat.smooth[-o,]
    gend.smooth.filt= gend.mat.smooth[-o,]
    tss.norm.smooth.filt = tss.norm.smooth[-o,]
    gend.norm.smooth.filt = gend.norm.smooth[-o,]
    
    tss.to.use = genes.large[-o]

    sample.outpath = file.path(outpath, sample.name)
        dir.create(sample.outpath, showWarnings=F)

    # ----------------------------------------------------- #
    # tss and gene start profile
    library(marray)
    cols.norm = maPalette(low="darkgray", mid="white", high="darkorange2", k=20)
    cols = colorRampPalette(c("lightgray", "darkblue"), bias=2)(50)
   
    a = rep(1:2, each=2)
    b=a
    a[seq(2, length(a), 2)] = 3
    indicator.matrix=cbind(b,a)

    DrawProfiles(mat.list = list(TSS=tss.mat.filt, Gene.End=gend.mat.filt, Random=random.mat), 
                 fact.list = list(NULL, NULL, NULL), 
                 indicator.matrix = indicator.matrix, 
                 name = paste(sample.name, "png", sep="."), 
                 outpath = sample.outpath, 
                 palette = c("darkorange","cornflowerblue","darkgray"), 
                 shift = c(2000, 20000), 
                 split = TRUE)
        
    # ----------------------------------------------------- #
    # function for plotting profiles with a given factor    
        PlotVariousProfiles = function(mat, mat.smooth, mat.norm, plot.outpath, factor, 
                                       genes, set.name, nset, shift){
 
            nfac = length(levels(factor))
            cat("Drawing profiles for :", set.name, "\n")   
           
            s = seq(1, nfac, nset)
            for(i in s){
                
                
                end = min(i+(nset-1),nfac)
                levs = levels(factor)[i:end, drop=F]
                print(levs)
                fac.ind = factor %in% levs
                n.levs=length(levs)
            
                tmp.fac = factor[fac.ind, drop=T]
				levels(tmp.fac) = levs
                mat.sel = mat[fac.ind,]
                mat.smooth.sel = mat.smooth[fac.ind,]
                mat.norm.sel = mat.norm[fac.ind,]
                genes.fac = genes[fac.ind]

                # ----------------------------- #
                # profiles
                cat("Drawing split profiles...\n")

                a = rep(1:n.levs, each=2)
                indicator.matrix=cbind(a,a) 
                indicator.matrix[seq(2, length(a), 2),2] = n.levs+1
                factor.palette = brewer.pal(n.levs, "Set1")[1:n.levs]
                DrawProfiles(mat.list=list(TSS=mat.sel, Random=random.mat), 
                             fact.list = list(tmp.fac, NULL), 
                             indicator.matrix=indicator.matrix, 
                             name=paste(set.name, i, "split","png", sep="."), 
                             outpath=plot.outpath, 
                             palette = c(factor.palette, "darkgray"), 
                             shift=shift, 
                             split=TRUE)
                
              
                DrawProfiles(mat.list=list(TSS=mat.sel, Random=random.mat), 
                             fact.list=list(tmp.fac, NULL), 
                             indicator.matrix=NULL, 
                             name=paste(set.name, i,"png", sep="."), 
                             outpath=plot.outpath, 
                             palette=c(factor.palette, "darkgray"), 
                             shift=shift, 
                             split=FALSE)
        }            
                # ----------------------------- #
                # heatmaps
                
                cat("Plotting the heatmaps...\n")
                width.ind = order(as.numeric(tmp.fac), width(genes.fac))
                DrawHeatmaps(mat=log10(mat.smooth.sel+1),
                             fact=tmp.fac, 
                             ord.fact=width.ind,
                             outpath=plot.outpath, 
                             name=paste(set.name, i, "Heat", "width", "png", sep="."), 
                             mat.cols=cols, 
                             key.cols=factor.palette)
                
                expr.ind = order(as.numeric(tmp.fac), values(genes.fac)$X0h.rpkm)
                DrawHeatmaps(mat=log10(mat.smooth.sel+1),
                             fact=tmp.fac, 
                             ord.fact=expr.ind,
                             outpath=plot.outpath, 
                             name=paste(set.name, i,"Heat", "expr", "png", sep="."), 
                             mat.cols=cols, 
                             key.cols=factor.palette)
                
                tss.smooth.filt.norm.2 = t(log10(apply( mat.norm.sel,1, function(x)(x-min(x))/(max(x)-min(x)))+1))
                width.ind = order(as.numeric(tmp.fac), width(genes.fac))
                DrawHeatmaps(mat=tss.smooth.filt.norm.2,
                             fact=tmp.fac, 
                             ord.fact=width.ind,
                             outpath=plot.outpath, 
                             name=paste(set.name, i,"Heat","norm", "width", "png", sep="."), 
                             mat.cols=cols.norm, 
                             key.cols=factor.palette)
            }

        }
    
          

               
    # -------------------------------------------------------- # 
    # bimodal
   
    bimodal = paste(values(tss.to.use)$H3k27me3_0h, sep="_", values(tss.to.use)$H3k4me3_0h)
    bimodal.f = as.factor(bimodal)
    levels(bimodal.f) = c("None","K4","K27","Bimodal")

    #  TSS
    bimodal.tss = file.path(sample.outpath, "Bimodal.TSS")
            dir.create(bimodal.tss, showWarnings=F)
    PlotVariousProfiles(mat = tss.mat.filt,
                        mat.smooth = tss.smooth.filt, 
                        mat.norm = tss.norm.smooth.filt,
                        plot.outpath = bimodal.tss, 
                        set.name = paste(sample.name,"TSS", "Bimodal", sep="."),
                        factor = bimodal.f,
                        genes=tss.to.use, 
                        nset=4,
                        shift=2000)

    
    # Gene.end
    bimodal.gend = file.path(sample.outpath, "Bimodal.GeneEnd")
            dir.create(bimodal.gend, showWarnings=F)
    PlotVariousProfiles(mat = gend.mat.filt,
                        mat.smooth = gend.smooth.filt, 
                        mat.norm = gend.norm.smooth.filt,
                        plot.outpath = bimodal.gend, 
                        set.name = paste(sample.name,"Gend", "Bimodal", sep="."),
                        factor = bimodal.f,
                        genes=tss.to.use,
                        nset=4,
                        shift=20000)


    
       


        
      
    # -------------------------------------------------------- # 
    # separate.biotypes

    # -------------------------------------------------------- # 
    # cpg island--histone modification(0h)
    bimodal.cpg = as.factor(paste(values(tss.to.use)$cpg.isl.prom, as.character(bimodal.f)))
    bimodal.tss.cpg = file.path(sample.outpath, "Bimodal.TSS.cpg")
        dir.create(bimodal.tss.cpg, showWarnings=F)
    PlotVariousProfiles(mat = tss.mat.filt,
                        mat.smooth = tss.smooth.filt, 
                        mat.norm = tss.norm.smooth.filt,
                        plot.outpath = bimodal.tss.cpg, 
                        set.name = paste(sample.name,"TSS", "Bimodal","cpg", sep="."),
                        factor = bimodal.cpg,
                        genes=tss.to.use, 
                        nset=4,
                        shift=2000)

    bimodal.gend.cpg = file.path(sample.outpath, "Bimodal.Gend.cpg")
        dir.create(bimodal.gend.cpg, showWarnings=F)
    PlotVariousProfiles(mat = gend.mat.filt,
                        mat.smooth = gend.smooth.filt, 
                        mat.norm = gend.norm.smooth.filt,
                        plot.outpath = bimodal.tss.cpg, 
                        set.name = paste(sample.name,"Gend", "Bimodal","cpg", sep="."),
                        factor = bimodal.cpg,
                        genes=tss.to.use, 
                        nset=4,
                        shift=20000)


        



    
    # -------------------------------------------------------- # 
    #TATA CpG　box--histone modification(0h)
    bimodal.tata = as.factor(paste(values(tss.to.use)$tata.box, as.character(bimodal.f)))
    bimodal.tss.tata = file.path(sample.outpath, "Bimodal.TSS.tata")
        dir.create(bimodal.tss.tata, showWarnings=F)
    PlotVariousProfiles(mat = tss.mat.filt,
                        mat.smooth = tss.smooth.filt, 
                        mat.norm = tss.norm.smooth.filt,
                        plot.outpath = bimodal.tss.tata, 
                        set.name = paste(sample.name,"TSS", "Bimodal","tata", sep="."),
                        factor = bimodal.tata,
                        genes=tss.to.use, 
                        nset=4,
                        shift=2000)

    bimodal.gend.tata = file.path(sample.outpath, "Bimodal.Gend.tata")
        dir.create(bimodal.gend.tata, showWarnings=F)
    PlotVariousProfiles(mat = gend.mat.filt,
                        mat.smooth = gend.smooth.filt, 
                        mat.norm = gend.norm.smooth.filt,
                        plot.outpath = bimodal.gend.tata, 
                        set.name = paste(sample.name,"Gend", "Bimodal","tata", sep="."),
                        factor = bimodal.tata,
                        genes=tss.to.use, 
                        nset=4,
                        shift=20000)



    # -------------------------------------------------------- # 
    # expression stratification 0h 72h
    expr.0h = log10(values(tss.to.use)$X0h.rpkm+1)
    expr.72h = log10(values(tss.to.use)$X72h.rpkm+1)

    expr.class.0h = paste(cut(expr.0h, breaks=c(0,1,4,5,6,max(expr.0h)), include.lowest=T), "0h",sep=".")
    expr.class.72h = paste(, cut(expr.72h, breaks=c(0,1,4,5,6,max(expr.72h)),	include.lowest=T), sep=".")
	expr.class = as.factor(paste(expr.class.0h, expr.class.72h, sep=' '))
	

    unique.ind = unique(expr.class.0h)
    expr.fact = as.factor(paste(expr.class.0h, expr.class.72h))
    table(expr.class.0h, expr.class.72h)

    expr.tss.outpath = file.path(sample.outpath, "Expression.0h.72h.TSS")
        dir.create(expr.tss.outpath, showWarnings=F)

    PlotVariousProfiles(mat = tss.mat.filt,
                        mat.smooth = tss.smooth.filt, 
                        mat.norm = tss.norm.smooth.filt,
                        plot.outpath =expr.tss.outpath, 
                        set.name = paste(sample.name,"TSS", "Expr", sep="."),
                        factor = expr.fact,
                        genes=tss.to.use, 
                        nset=4,
                        shift=2000)



    expr.gend.outpath = file.path(sample.outpath, "Expression.0h.72h.Gend")
        dir.create(expr.gend.outpath, showWarnings=F)

    PlotVariousProfiles(mat = gend.mat.filt,
                    mat.smooth = gend.smooth.filt, 
                    mat.norm = gend.norm.smooth.filt,
                    plot.outpath =expr.gend.outpath, 
                    set.name = paste(sample.name,"Gend", "Expr", sep="."),
                    factor = expr.fact,
                    genes=tss.to.use, 
                    nset=4,
                    shift=20000)

    # -------------------------------------------------------- # 
    # bimodal 0h - bimodal 72h
    bimodal.0h = as.factor(paste(values(tss.to.use)$H3k27me3_0h, sep="_", values(tss.to.use)$H3k4me3_0h))
	levels(bimodal.0h) = paste(c("No","K4","K27","Bi"),"0h",  sep='.')
    bimodal.72h = as.factor(paste(values(tss.to.use)$H3k27me3_72h, sep="_", values(tss.to.use)$H3k4me3_72h))
	levels(bimodal.72h) = paste(c("No","K4","K27","Bi"), "72h", sep='.')
	levs = apply(expand.grid(levels(bimodal.0h), levels(bimodal.72h)), 1, paste, collapse=" ")

    
    bimodal.fact = factor(paste(as.character(bimodal.0h), as.character(bimodal.72h)))
	levels(bimodal.fact) = c(levels(bimodal.fact), setdiff(levs, levels(bimodal.fact)))

    bimodal.tss.outpath = file.path(sample.outpath, "Bimodal.0h.72h.TSS")
        dir.create(bimodal.tss.outpath, showWarnings=F)
    
    PlotVariousProfiles(mat = tss.mat.filt,
                    mat.smooth = tss.smooth.filt, 
                    mat.norm = tss.norm.smooth.filt,
                    plot.outpath =bimodal.tss.outpath, 
                    set.name = paste(sample.name,"TSS", "BimodCond", sep="."),
                    factor = bimodal.fact,
                    genes=tss.to.use, 
                    nset=4,
                    shift=2000)


    bimodal.gend.outpath = file.path(sample.outpath, "Bimodal.0h.72h.Gend")
        dir.create(bimodal.gend.outpath, showWarnings=F)

    PlotVariousProfiles(mat = gend.mat.filt,
                    mat.smooth = gend.smooth.filt, 
                    mat.norm = gend.norm.smooth.filt,
                    plot.outpath = bimodal.gend.outpath, 
                    set.name = paste(sample.name,"Gend", "BimodCond", sep="."),
                    factor = bimodal.fact,
                    genes=tss.to.use, 
                    nset=4,
                    shift=20000)

}
