ScalerLarge = function(a, len){
    
    s = unique(ceiling(seq(1, length(a)-1, length.out=len)))
    win = floor(length(a)/len)
    win.starts = ceiling(s-win/2)
    win.starts[1] = 1
    win.ends = ceiling(s+win/2)
    win.ends[length(win.ends)] = length(a)
    v = round(viewMeans(Views(a, win.starts, win.ends)))
    return(v)
}

library(IRanges)

v = c(1, 1, 1, 2, 2, 2, 1, 1, 1, 1, 0, 0, 0, 0, 3, 3, 3, 4, 4, 4)
a = Rle(v)

s = unique(ceiling(seq(1, length(a)-1, length.out=len)))
    win = floor(length(a)/len)
    win.starts = ceiling(s-win/2)
    win.starts[1] = 1
    win.ends = ceiling(s+win/2)
    win.ends[length(win.ends)] = length(a)
    v = round(viewMeans(Views(a, win.starts, win.ends)))

outpath= "/common/SHARED/vfranke/Fugaku_ChipSeq/Data"
library(BSgenome.Mmusculus.UCSC.mm9)
library(stringr)
seqlen = seqlengths(Mmusculus)
seqlen = seqlen[!str_detect(names(seqlen), "random")]
wins = MakeTillingWindow(seqlengths(Mmusculus), 200)
wins = wins[!str_detect(as.vector(seqnames(wins)), "random")]
wins.s = split(wins, seqnames(wins))

chrs = as.character(unique(seqnames(wins)))
outfile = file.path(outpath, "gc.mm9.bedGraph")
if(file.exists(outfile))
    file.remove(outfile)
for(i in chrs){
    cat(i, "\n")
    wins.tmp = wins.s[[i]]
    v = Views(Mmusculus[[i]],start(wins.tmp), end(wins.tmp))
    gc = oligonucleotideFrequency(v, width=2)
    d = data.frame(i, start(wins.tmp), end(wins.tmp), gc[,10])
    write.table(d, outfile, col.names=F, row.names=F, quote=F, sep="\t", append=T)
}

