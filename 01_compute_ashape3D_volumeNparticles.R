options(rgl.printRglwidget = TRUE,scipen=9999)
library(alphashape3d)
library(rgl)
args = commandArgs(trailingOnly=TRUE)

# Values of alpha
factor = 0.5 * sqrt(2)
alphas = seq(factor, factor*(as.integer(20/factor))+1, factor)
#alphas = seq(1, as.integer(20/factor)+1, 1.0)
print(alphas)
#quit()
triplets <- list(c("A","B","C"),c("B","A","C"),c("C","A","B"))

x <- read.table("_xyzAall")
x <- as.matrix(x)
x[] <- as.double(x)

y <- read.table("_xyzBall")
y <- as.matrix(y)
y[] <- as.double(y)

z <- read.table("_xyzCall")
z <- as.matrix(z)
z[] <- as.double(z)

for(triplet in triplets)
{
    print(paste(triplet,collapse=" "))
    if(triplet[1]=="A")
    {
        x1 <- x
	x2 <- y
	x3 <- z
    }
    if(triplet[1]=="B")
    {
        x1 <- y
	x2 <- x
	x3 <- z
    }
    if(triplet[1]=="C")
    {
        x1 <- z
	x2 <- x
	x3 <- y
    }

    # 3D alpha-shape for all alphas for chain A
    print(x1)
    ashape3d.obj <- ashape3d(x1, alpha = factor*alphas, pert=TRUE)
    #General Position. Throughout this article we assume that the points of S
    #are in general position. For the time being, this means that no 4 points lie on
    #a common plane; no 5 points lie on a common sphere; and for any fixed (Y, the
    #smallest sphere through any 2,3, or 4 points of S has a radius different from
    #LY. The general-position assumption will later be extended when convenient
    #(see Section 5.3). It simplifies forthcoming definitions, discussions, and algorithms and is justified by a programming technique known as SOS
    #[Edelsbrunner and Miicke 19901. This method simulates an infinitesimal
    #perturbation of the points on the level of geometric predicates and relieves
    #the programmer from the otherwise necessary case analysis (see Section 6.2).
    #eps = 1e-09
    outFile = paste0(triplet[1],"volume.tab")
    volumes = volume_ashape3d(ashape3d.obj, byComponents = FALSE, indexAlpha = 1:length(alphas))
    df <- data.frame(alphas,volumes)
    write.table(df,file=outFile,row.names = F,quote = FALSE, col.names=F, sep="\t", append=T)

    for(j in c(2,3))
    {
        outFile = paste0(triplet[j],"particles_in_",triplet[1],"_list.tab")
	outFile1 = paste0(triplet[j],"particles_in_",triplet[1],".tab")	
    	for(i in seq(1,length(alphas),1))
    	{
            #print(outFile1)
            
	    if(j==2)
	    {
    	        inside = as.list(x2[inashape3d(ashape3d.obj, indexAlpha = i, x2),])
            } else {
                inside = as.list(x3[inashape3d(ashape3d.obj, indexAlpha = i, x3),])
            }
	    #print(head(inside))
    	    #print(alphas[[i]])
    	    #print(typeof(inside))
    	    #print(length(inside))
    	    df <- data.frame(alphas[[i]],length(inside)/3)
    	    write.table(df,file=outFile1,row.names = F,quote = FALSE, col.names=F, sep="\t", append=T)

    	    next
    	    if(dim(inside)[[1]] == 0){next;}
    	    df <- data.frame(alphas[[i]],inside)
    	    colnames(df) <- c("replica","timestep","alpha","x","y","z")
    	    #print(head(df))    
    	    write.table(df,file=outFile,row.names = F,quote = FALSE, col.names=F, sep="\t", append=T)
        }
    }
}
quit()
