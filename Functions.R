# INFO: Holds basic functions for the ChipSeqProject
# DATE: 28.11.2011.
# AUTH: ABC


# -------------------------------------------------------------- #
# {1}
    ### loads an RData file and assigns it to the name variable
    Assigner = function(`_path`, `_name`){
        if(! is.character(`_path`) | !is.character(`_name`))
        stop('Both arguments should be characters!')
        load(`_path`)
        assign(`_name`, get(ls()[1]), parent.frame())
    }


# -------------------------------------------------------------- #
# {2}
    # Extends a region around a given genomic mark (e.g. TSS)
    # As input takes a data frame - first three columns have to bee
    # chr   start end
    # strand is optional
    RegionExtender = function(bed, downstream = 1000, upstream = 1000, strand = F){
    
        if(strand == T){
        
            ### checks if the strand column exists
            if (length(grep("strand", names(bed))) == 0){
                cat('Please specify the strand column!\n')
                return(NULL)
            }
        
            bed.m = bed$strand == '-'
            bed[bed.m,2] = bed[bed.m,2] - downstream
            bed[bed.m,3] = bed[bed.m,3] + upstream
        
            bed.p = bed$strand == '+'
            bed[bed.p,2] = bed[bed.p,2] - upstream
            bed[bed.p,3] = bed[bed.p,3] + downstream
            colnames(bed)[2:3] = c('start','end')
            return(bed)
        }
        else if(strand == F){
            bed[,2] = bed[,2] - upstream
            bed[,3] = bed[,3] + downstream
            colnames(bed)[2:3] = c('start','end')
            return(bed)
        }
    }


# -------------------------------------------------------------- #
# {3}
    ### takes a bed file and designates start and end as same coordinates, based on the strand
    StartStrandMaker = function(bed){
    
        ### for + strand, end = start
        strand = grep('strand', names(bed))
        bed[bed[, strand] == '+', 3] = bed[bed[,strand] == '+', 2]
        ### for - strand, start = end
        bed[bed[,strand] == '-', 2] = bed[bed[,strand] == '-', 3]
        return(bed)
    }


# -------------------------------------------------------------- #
# {4}
    # Converts a bed formatted data frame to a GRanges object
    BedToGRanges = function(bed, values=FALSE, seqlen=NULL){
            
        library(stringr)
        library(GenomicRanges)
        # sorts the data frame by chromosome - position
        bed = bed[order(bed[,1], bed[,2]),]
        
        # checks whether the bed file contains the strand column
        if(any(str_detect(names(bed), 'strand'))){
            strand=as.character(bed[,str_detect(names(bed), 'strand')])
        }else{  
            strand='*'
        }
        
        # checks whether you have provided the chromosome lengths
        if(!is.null(seqlen)){
            bed = RegionCorrector(bed, seqlen)
            ranges = GRanges(
                             seqnames=as.character(bed[,1]),
                             ranges=IRanges(start=as.numeric(bed[,2]), end=as.numeric(bed[,3])),
                             strand=strand,
                             seqlengths=seqlen
                            )
        }else{  
             ranges = GRanges(
                              seqnames=as.character(bed[,1]),
                              ranges=IRanges(start=as.numeric(bed[,2]), end=as.numeric(bed[,3])),
                              strand=strand
                              )
        }
        
        # fills the data frame with annotation columns
        if(values == TRUE){
            col.ids = setdiff(names(bed), c( "seqnames", "ranges", "strand", "seqlengths", "start", "end", "width",  "element", "chr"))
            values(ranges) = bed[col.ids]
        }
            
        return(ranges)
    }

# -------------------------------------------------------------- #
#{5}
    # Make tilling window over the genome
    MakeTillingWindow = function(seqlen, window.size = 2000){
        
        tile.list = list()
        for(i in names(seqlen)){
            cat(i,"\r")
            start = seq(1, seqlen[i], by=window.size)
            g = GRanges(i, IRanges(start=start, width=window.size), strand="+")
            g = g[end(g) < seqlen[i]]
            seqlevels(g) = names(seqlen)
            seqlengths(g) = seqlen
            tile.list[[i]] = g
        
        }
        cat("\n")
        tiles = Reduce(c,tile.list)
        return(tiles)
    }

# -------------------------------------------------------------- #
#{6}
    # Creates random regions and takes only those that do not overlap any segment
    CreateRandomRegions = function(seqlen, window.size, ranges=NULL, num.of.rand.reg=5000, n=2, k=1){
        
        wins = MakeTillingWindow(seqlen, (window.size*n+k))
        if(!is.null(ranges)){
            cat("Removing the random regions that overlap designated regions...\n")
            overlap = as.matrix(findOverlaps(wins, ranges))
            wins = wins[-overlap[,1]]
        }
        rand.reg = wins[sample(1:length(wins), num.of.rand.reg)]
        rand.reg.s = split(rand.reg, seqnames(rand.reg))
        return(rand.reg.s)

    }   

# -------------------------------------------------------------- #
#{7}
    #Takes a granges object and sets the coordinate to start/end based on the strand
    MakeViewPoint = function(granges, viewpoint="start"){
        
        if(!class(granges) == "GRanges")
            stop("Input must be a valid GRanges object")

        if(!viewpoint %in% c("start","end") | length(viewpoint) >1)
            stop("viewpoint argument can only take values of start or end")
        
        amb.strand.ind = as.vector(strand(granges) == "*")
        if(any(amb.strand.ind)){
            cat("Setting ambiguous strand to +\n")
            strand(granges)[amb.strand.ind] = "+"
        }
            
        pos.strand = as.vector(strand(granges) == "+")
        if(viewpoint == "start"){
            cat("Setting the viewpoint to start...\n")
            end(granges)[pos.strand]    = start(granges)[pos.strand]
            start(granges)[!pos.strand] = end(granges)[!pos.strand] 
        }else if(viewpoint == "end"){
            cat("Setting the viewpoint to end...\n")
            end(granges)[!pos.strand]  = start(granges)[!pos.strand]
            start(granges)[pos.strand] = end(granges)[pos.strand] 
        }

        return(granges)
    }

# -------------------------------------------------------------- #
#{8}
    #Expands a granges object around the viewpoint
    ExpandGranges = function(granges, upstream=2000, downstream=2000, strand = F){

        if(!class(granges) == "GRanges")
            stop("Input must be a valid GRanges object")

        if(!is.numeric(upstream))
            stop("Upstream must be numeric")

        if(!is.numeric(downstream))
            stop("Downstream must be numeric")

        if(!is.logical(strand))
            stop("strand must be logical")

                
        if(strand == FALSE){
            cat("Expanding the unstranded regions...\n")
            start(granges) = start(granges) - upstream
            end(granges) = end(granges) + downstream
        }

        if(strand == TRUE){
            cat("Expanding the stranded regions...\n")
            strand.ind = as.vector(strand(granges) == "+")
            start(granges)[strand.ind] = start(granges)[strand.ind] - upstream
            end(granges)[strand.ind] = end(granges)[strand.ind] + downstream
            
            start(granges)[!strand.ind] = start(granges)[!strand.ind] - downstream
            end(granges)[!strand.ind] = end(granges)[!strand.ind] + upstream
        }
        return(granges)

    }


# -------------------------------------------------------------- #
# {9}
    # fives the index of granges instances which fell of the chromosomes
    OffChromosomeRegions = function(granges, seqlen){

        ind = rep(FALSE,length(granges))
        for(i in names(seqlen)){
            cat(i,"\r")
            ind[as.vector(seqnames(granges) == i) & 
                as.vector(end(granges)) < seqlen[i] & 
                start(granges) > 0] = TRUE
        }
        return(ind)
    }

