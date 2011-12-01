# INFO: Functions for getting profile infromation from coverage vectors
# DATE: 01.12.2011.
# AUTH: ABC

# -------------------------------------------------------------- #
# {1}
#
Coverage2Profiles = function(cov, ranges, len=1000){
    
    list.mat = GetProfiles(cov, ranges)
    
    # if all elements of the list are matrices of the same size then it concatenates them into a data frame and reverses the profile on the - strand
    # if all elements of the list are lists (views from different sizes) - it does the smoothing, to get them into same size matrices, concatenates them and converts them into a data.frame
    if(all(sapply(list.mat, class)) == "matrix"){
        if(!nrow(unique(t(sapply(list.mat, dim)))) == 1)
            stop("All matrices do not have the same dimension")
        mat = MatListToDataFrame(list.mat, ranges)
    
        
    }else if(all(sapply(list.mat, class)) == "list"){
        if(any(sapply(list.mat[[i]], length) < len))
        stop("All views have to have the length greater than the smoothing factor")
        mat = ListListToDataFrame(list.mat, ranges, len=len)
        
        
    }else{
        stop("Not all elements of the list correspond to the same class")
    }
    
    return(mat)
    
}



# -------------------------------------------------------------- #
# {2}
#gets the profile from the coverage file and the GRanges object
GetProfiles = function(cov, ranges){
    
    if(!class(cov) == "SimpleRleList")
        stop("Coverage needs to be a SimpleRleList object")
    if(!class(ranges) == "GRangesList")
        stop("Ranges needs to be a granges object")
    
    cov = cov[order(names(cov))]
    ranges = ranges[order(names(ranges))]  
    
    if(! all(names(sample.coverage) == names(windows.s)))
        stop("List names do not match")
    
    
    # getts the views around the tss
    cat("Getting the views...\n")
    v = Views(sample.coverage, ranges(windows.s))
    vm = lapply(v, function(x)viewApply(x, as.vector))
    
    return(vm)
}


# -------------------------------------------------------------- #
# {3}
    # Converts a list of matrices to a data frame and reverses the profile on the - strand
    MatListToDataFrame = function(mat.list, ranges){
    
        cat("Converting the list to the matrix...\n")
        mat = data.frame(do.call(rbind, lapply(vm, t)))
        # reverses the profiles on the negative strand
        cat("Reversing the matrices...\n")
        strand.ind = unlist(lapply(ranges, strand)) == "-"
        if(nrow(mat) != length(strand.ind))
            stop("Strand and matrix sizes do not match")
        mat[strand.ind,] = t(apply(mat[strand.ind,], 1, rev))
        return(mat)
    }



# -------------------------------------------------------------- #
# {4}
    # Converts a list of lists to a data frame and reverses the profiles on the - strand
    ListListToDataFrame = function(list.mat, ranges, len){

        # gets the strand information from the ranges
        strand.ind = lapply(ranges, function(x)strand(x) == "-"))
        # reverses each element that is on the - strand
        cat("Reversing the lists...\n")
        for(i in 1:length(list.mat)){
            list.mat[[i]][strand.ind[[i]]] = lapply(list.mat[[i]][strand.ind[[i]]], rev)
        }
        
        cat("Smoothing the profiles...\n")
        for(i in 1:length(list.mat)){
            if(any(sapply(list.mat[[i]], length) < len))
                stop("All views have to have the length greater than the smoothing factor")
            list.mat[[i]] = lapply(list.mat[[i]], ScalerLarge, len)
        }
        mat = data.frame(do.call(rbind, list.mat))
        return(mat)
        
    }




# -------------------------------------------------------------- #
# {5}
    # takes an Rle object and smooths it
    ScalerLarge = function(a, len){
        
        s = unique(ceiling(seq(1, length(a), length.out=len)))
        win = floor(length(a)/len - 1)
        win.starts = ceiling(s-win/2)
        win.starts[1] = 1
        win.ends = ceiling(s+win/2)
        win.ends[length(win.ends)] = length(a)
        v = round(viewMeans(Views(a, win.starts, win.ends)))
        return(v)
    }
