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




