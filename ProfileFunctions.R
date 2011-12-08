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
    
        
    }else if(all(sapply(mat.list, class) == "list")){
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
    mind = sapply(vm, is.matrix)
    if(!all(mind)){
        vm[!mind] = lapply(vm[!mind], function(x)matrix(x, ncol=1))
    }
    
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
    ListListToDataFrame = function(mat.list, ranges, len){

        # gets the strand information from the ranges
        strand = strand(ranges)
        strand = lapply(strand, as.vector)
        strand = lapply(strand, function(x)x[x == "*"] = "+")
        strand.ind = lapply(strand, function(x)x == "-")
        # reverses each element that is on the - strand
        cat("Reversing the lists...\n")
        for(i in 1:length(mat.list)){
            if(any(strand.ind[[i]]))
                mat.list[[i]][strand.ind[[i]]] = lapply(mat.list[[i]][strand.ind[[i]]], rev)
        }
        
        cat("Smoothing the profiles...\n")
        for(i in 1:length(mat.list)){
            cat(i,"\r")
            if(any(sapply(mat.list[[i]], length) < len))
                stop("All views have to have the length greater than the smoothing factor")
            mat.list[[i]] = lapply(mat.list[[i]], ScalerLarge, len)
        }
        mat = do.call(rbind, lapply(mat.list, function(x)data.frame(do.call(rbind, x))))
        return(mat)
        
    }




# -------------------------------------------------------------- #
# {5}
    # takes an Rle object and smooths it
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


# -------------------------------------------------------------- #
# {6}
    #Takes a matrix and a factor and draws the profiles based on the factor
    DrawProfiles = function(mat.list, fact.list=NULL, indicator.matrix=NULL, name, outpath, palette, shift, split=FALSE){
        
        if(is.null(mat.list))
            stop("The list of matrices need to be designated")
        if(is.null(name))
            stop("The output name needs to be designated")
        if(is.null(outpath))
            stop("The output directory needs to be designated")
        
        profile.list = CalculateColMeans(mat.list, fact.list)
        ProfilePlotter(l=profile.list, 
                       m=indicator.matrix, 
                       split=split, 
                       outpath=outpath, 
                       name=name, 
                       palette=palette, 
                       shift=shift)
            
    }

# -------------------------------------------------------------- #
# {7}
    # Takes a list of matrices and a list of factors and returns a list of caluculated colmeans for each matrix for each factor
    CalculateColMeans = function(mat.list, fact.list=NULL){
        

        if(any(sapply(names(mat.list), is.null)))
            stop("All mat.list elements have to have designated names")
        
        cat("Calculating cumulative profiles\n")
        
        colmeans.list = list()
        for(i in 1:length(mat.list)){
                
            cat(i,"\r")
            if(!is.null(fact.list) & !is.null(fact.list[[i]])){
                if(nrow(mat.list[[i]]) != length(fact.list[[i]]))
                    stop("Each profile does not have a designated factor")
                colmeans.tmp = lapply(split(mat.list[[i]], fact.list[[i]]), colMeans)
                colmeans.list = c(colmeans.list, colmeans.tmp)
                
            }else{
                colmeans.list = c(colmeans.list, list(colMeans(mat.list[[i]])))
                names(colmeans.list)[length(colmeans.list)] = names(mat.list)[i]
                
            }
        }
        return(colmeans.list)
    }


# -------------------------------------------------------------- #
# {8}
    # Takes a list of vectors (mean coverage over a window), and a designator matrix which tells which profiles to plot on the same plot
        ProfilePlotter = function(l = NULL, m = NULL, split=FALSE, outpath, name, palette, shift=0){
                
            if(!is.list(l))
                stop("l needs to be a list")        
                
                
            if(split == FALSE)
                m = cbind(1,1:length(l))
            if(split == TRUE & is.null(m))
                stop("When split is true indicator matrix must be given")
                
            if(!is.matrix(m))
                stop("m needs to be matrix")
            if(any(!m %in% 1:length(l)))
                stop("matrix has invalid combinations for plotting")
            nclass = unique(m[,1])
            height = max(c(1000, 350*(length(nclass))))
            png(file.path(outpath, name), width=1000, height=height)
                
            if(split == TRUE){
                par(mfrow=c(length(nclass), 1),cex=max(1.25, 0.30*length(nclass)))
                cat("Drawing the profiles with split...\n")
            }
            if(length(shift) != length(nclass))
                shift = rep(shift, length(nclass))
        
            for(i in nclass){
                        
                m.tmp = m[m[,1] == i,]
                l.tmp = l[m.tmp[,2]]
                palette.tmp = palette[m.tmp[,2]]
                x = 0:(length(l[[i]])-1) - shift[i]
                for(j in 1:length(l.tmp)){
                    if(j == 1){
                        ylim = c(min(unlist(l.tmp, use.names=F)), max(unlist(l.tmp, use.names=F)))
                        plot(x, l.tmp[[j]], col=palette.tmp[j], type="l", lwd=2, ylab="mean.coverage", xlab="position", main=names(l.tmp)[j], ylim=ylim)
                    }else{
                        lines(x, l.tmp[[j]], col=palette.tmp[j], lwd=2)
                    }
                }
                legend("topright", fill=palette.tmp, legend=names(l.tmp))
            }
            dev.off()
        }
