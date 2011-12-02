# INFO: Functions for getting profile infromation from coverage vectors
# DATE: 01.12.2011.
# AUTH: ABC

# -------------------------------------------------------------- #
# {1}
#
Coverage2Profiles = function(cov, ranges, len=1000){
    
    if(!class(cov) == "SimpleRleList")
        stop("Coverage needs to be a SimpleRleList object")
    if(!class(ranges) == "GRangesList")
        stop("Ranges needs to be a granges object")
    
    cov = cov[order(names(cov))]
    ranges = ranges[order(names(ranges))]  
    
    if(! all(names(cov) == names(ranges)))
        stop("List names do not match")
    
    mat.list = GetProfiles(cov, ranges)
    mat.list = mat.list[!is.null(sapply(mat.list, nrow))]
    
    # if all elements of the list are matrices of the same size then it concatenates them into a data frame and reverses the profile on the - strand
    # if all elements of the list are lists (views from different sizes) - it does the smoothing, to get them into same size matrices, concatenates them and converts them into a data.frame
    if(all(sapply(mat.list, class) == "matrix")){
         mat = MatListToDataFrame(mat.list, ranges)
    
        
    }else if(all(sapply(mat.list, class)) == "list"){
        if(any(sapply(mat.list[[i]], length) < len))
        stop("All views have to have the length greater than the smoothing factor")
        mat = ListListToDataFrame(mat.list, ranges, len=len)
        
        
    }else{
        stop("Not all elements of the list correspond to the same class")
    }
    
    return(mat)
    
}



# -------------------------------------------------------------- #
# {2}
#gets the profile from the coverage file and the GRanges object
GetProfiles = function(cov, ranges){
    
    # getts the views around the tss
    cat("Getting the views...\n")
    v = Views(cov, ranges(ranges))
    vm = lapply(v, function(x)viewApply(x, as.vector))
    
    return(vm)
}


# -------------------------------------------------------------- #
# {3}
    # Converts a list of matrices to a data frame and reverses the profile on the - strand
    MatListToDataFrame = function(mat.list, ranges){
    
        cat("Converting the list to the matrix...\n")
        mat = data.frame(do.call(rbind, lapply(mat.list, t)))
        # reverses the profiles on the negative strand
        cat("Reversing the matrices...\n")
        strand = unlist(lapply(ranges, function(x)as.vector(strand(x))))
        strand[strand == "*"] = "+"
        strand.ind = strand == "-"
        if(nrow(mat) != length(strand.ind))
            stop("Strand and matrix sizes do not match")
    
        if(any(strand.ind))
            mat[strand.ind,] = t(apply(mat[strand.ind,], 1, rev))
        return(mat)
    }



# -------------------------------------------------------------- #
# {4}
    # Converts a list of lists to a data frame and reverses the profiles on the - strand
    ListListToDataFrame = function(list.mat, ranges, len){

        # gets the strand information from the ranges
        strand = strand(ranges)
        strand = lapply(strand, function(x)x[x == "*"] = "+")
        strand.ind = lapply(strand, function(x)x == "-")
        # reverses each element that is on the - strand
        cat("Reversing the lists...\n")
        for(i in 1:length(list.mat)){
            if(any(strand.ind[[i]]))
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


# -------------------------------------------------------------- #
# {6}
    #Takes a matrix and a factor and draws the profiles based on the factor
    DrawProfiles = function(mat.list, fact.list=NULL, indicator.mat=NULL name, outpath, palette, shift, split=FALSE){
        
        if(is.null(mat.list))
            stop("The list of matrices need to be designated")
        if(is.null(name))
            stop("The output name needs to be designated")
        if(is.null(outpath))
            stop("The output directory needs to be designated")
        profile.list = CalculateColMeans(mat.list, fact.list)
        ProfilePlotter(l=profile.list, m=indicator.matrix, split=split, outpath=outpath, name=name, cols=cols, shift=shift))
            
    }

# -------------------------------------------------------------- #
# {7}
    # Takes a list of matrices and a list of factors and returns a list of caluculated colmeans for each matrix for each factor
    CalculateColMeans = function(mat.list, fact.list=NULL){

        cat("Calculating cumulative profiles\n")
        colmeans.list = list()
        for(i in 1:length(mat.list)){
                
            cat(i/length(mat.list),"\r")
            list.len = length(colsum.list)
            if(!is.null(fact.list) && !is.null(fact.list[[i]])){
                colmeans.tmp = lapply(split(mat.list[[i]], fact.list[[i]]), colMeans)
                colmeans.list = c(colmeans.list, lapply(split(mat.list[[i]], fact.list[[i]]), colMeans))
                
            }else{
                colmeans.list = c(colmeans.list, (mat.list[[i]]))
                
            }
        }
        return(colmeans.list)
    }


# -------------------------------------------------------------- #
# {8}
    # Takes a list of vectors (mean coverage over a window), and a designator matrix which tells which profiles to plot on the same plot
    ProfilePlotter = function(l = NULL, m = NULL, split=F, outpath, name, cols, shift=0){

        if(!is.list(l))
            stop("l needs to be a list")        
        
        
        if(split == FALSE){
                
            cat("Drawing the profiles without the split...\n")
            png(file.path(outpath, name), width=1200, height=800)
            for(i in 1:length(l))
                
                if(i == 1){
                    ylim = c(min(unlist(l, use.names=F)), max(unlist(l, use.names=F)))
                    x = 1:length(l[[i]]) - shift
                    plot(x, l[[i]], cols[[i]], type="l", lwd=2, ylab="mean.coverage", xlab="position", main=name, ylim=ylim)
                }else{
                    lines(x, l.tmp[[j]], col=cols[[i]], lwd=2)
                }
            }
            legend("topright", fill=cols, legend=names(l))
            
        
        }else{

            if(!is.matrix(m))
                stop("m needs to be matrix")
            if(any(!m %in% 1:length(l)))
                stop("matrix has invalid combinations for plotting")
            cat("Drawing the profiles with split...\n")
            nclass = unique(m[,1])
            png(file.path(outpath, name), width=1000, height=400*nclass)
            
            par(mfrow=c(nclass, 1)
            for(i in 1:nclass){
                
                m.tmp = m[m[,1] == i,]
                l.tmp = l[m.tmp[,2]]
                for(j in m.tmp[,2]){
                    if(j == 1){
                        ylim = c(min(unlist(l.tmp, use.names=F)), max(unlist(l.tmp, use.names=F)))
                        x = 1:length(l[[i]]) - shift
                        plot(x, l.tmp[[j]], cols[[j]], type="l", lwd=2, ylab="mean.coverage", xlab="position", main=names(l.tmp)[j], ylim=ylim)
                    }else{
                        lines(x, l.tmp[[j]], col=cols[[j]], lwd=2)
                }
            }
        }
        dev.off()
    }