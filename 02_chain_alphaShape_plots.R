options(scipen=9999)
library(dplyr)
library(ggplot2)
library(ggpubr)
library("ggrepel")
library(scales)
library("RColorBrewer")
pi=3.14159265359

#Radius of confinement = 10 Lattice units 
radius=10.0 # Radius of the confining sphere in sigma
totVolume=(4./3.*pi*radius*radius*radius)
print(paste0("Total nuclear volume ",totVolume))

factor = 0.5 * sqrt(2)
alphas = seq(factor, factor*(as.integer(20/factor))+1, factor)
#alphas = seq(1, as.integer(20/factor)+1, 1.0)
print(alphas)
#

totBeads=2000
print(paste0("Total beads ",totBeads))

data <- read.table("chain_alphaShape_volume.tab", header=F)
data <- data[c(1,2,3)]
colnames(data) <- c("alpha","timestep","alphaShapeVol")
data$alphaShapeVol <- data$alphaShapeVol / totVolume
print(head(data))

colors <- rainbow(length(unique(data$alpha)))
#colors <- scale_colour_brewer(palette="Reds")(length(unique(data$alpha)))
print(paste0("Number of alphas ",length(unique(data$alpha))))

print(paste0("alphaShape"," ","volume"))
tolerance=0.005
outFileD <- "alphaShapeVol_vs_time.tab"
string <- data.frame("timestep","alphaPlateau","alphaShapeVolPlateau")
write.table(string,file=outFileD,row.names = F,quote = FALSE, col.names=F, sep="\t", append=T)

for(time in seq(0,100,1))
{
    # Plot the ashapeVol per replica for some values of timestep
    outFile <- paste0("chain_alphaShapeVol_vs_alpha_time_",time,".pdf")
    if(!file.exists(outFile))
    {
        print(paste0("outFile ",outFile))
	dataPlot <- data[data$timestep == time & data$alpha <= 15,]
	print(head(dataPlot))

	p <- ggplot(dataPlot, aes(x=factor(alpha), y=alphaShapeVol, fill="red")) +
     	     geom_violin( position = position_dodge(1), trim=T, scale = "width") +
	     #scale_fill_manual(values=colors) +	     
	     geom_boxplot(position = position_dodge(1), fill="white", color="black", width=0.1, outlier.shape = NA) +     
	     ylim(0.0,1.0) +
	     theme(panel.background = element_rect(fill = NA),
                 panel.grid.major.x = element_blank(),
                 panel.grid.major.y = element_blank(),
            	 axis.title = element_text(face="bold",size=24),
		 axis.text = element_text(face="bold",size=12),
		 axis.text.x = element_text(angle=60, hjust=1),
            legend.position = "none") +
  	    labs(x="alpha",y="Fract. alpha-shape-volume",fill="") +     
	    guides(fill="none") #+
	    #geom_signif(comparisons = my_comparisons, y_position=c(ymax*(1+0.05),ymax*(1+0.05))) +
	    #geom_text(data=summ, aes(y=score, label = as.character(n)), color="black")

    	pdf(outFile)
    	print(p)
    	dev.off()
    }
}



for(alpha in alphas)
{
    # Plot the ashapeVol per replica for some values of timestep
    #outFile <- paste0("chain_alphaShapeVol_vs_time_alpha_",alpha,".pdf")
    outFile <- paste0("chain_alphaShapeVol_vs_time_tMax_20_alpha_",alpha,".pdf")    
    if(!file.exists(outFile))
    {
        print(paste0("outFile ",outFile))
	print(paste0("alpha ",alpha))
	#dataPlot <- data[-0.1 < (data$alpha - alpha) & (data$alpha - alpha) < 0.1,]
	dataPlot <- data[-0.1 < (data$alpha - alpha) & (data$alpha - alpha) < 0.1 & data$timestep<=20,]	
	print(head(data))
	print(head(dataPlot))
	#quit()
	
	p <- ggplot(dataPlot, aes(x=factor(timestep), y=alphaShapeVol, fill="red")) +
     	     geom_violin( position = position_dodge(1), trim=T, scale = "width") +
	     #scale_fill_manual(values=colors) +	     
	     geom_boxplot(position = position_dodge(1), fill="white", color="black", width=0.1, outlier.shape = NA) +     
	     #ylim(0.0,1.0) +
	     theme(panel.background = element_rect(fill = NA),
                 panel.grid.major.x = element_blank(),
                 panel.grid.major.y = element_blank(),
            	 axis.title = element_text(face="bold",size=24),
		 axis.text = element_text(face="bold",size=12),
		 axis.text.x = element_text(angle=60, hjust=1),
            legend.position = "none") +
  	    labs(x="time",y="Fract. alpha-shape-volume",fill="") +     
	    guides(fill="none") #+
	    #geom_signif(comparisons = my_comparisons, y_position=c(ymax*(1+0.05),ymax*(1+0.05))) +
	    #geom_text(data=summ, aes(y=score, label = as.character(n)), color="black")

    	#p <- ggplot() +
	#    geom_violin(dataPlot, mapping=aes(x=factor(timestep), y=alphaShapeVol, color=rainbow(3)[[1]], alpha=0.2)) +	 
        #    ylim(0.0,1.0) +
        #    theme(panel.background = element_rect(fill = NA),
        #         panel.grid.major.x = element_blank(),
        #         panel.grid.major.y = element_blank(),
        #    	 axis.title = element_text(face="bold",size=24),
	#	 axis.text = element_text(face="bold",size=18),
        #    legend.position = "none") +
            #labs(x="time",y="Fract. alpha-shape-volume",fill="") +
        #    labs(x="time",y="Alpha-shape-volume",fill="") +     	    
            #geom_text(data=stats,mapping=aes(x=alpha,y=0, label = n), color="black") +
        #    guides(fill="none")

    	pdf(outFile)
    	print(p)
    	dev.off()
    }
}

print(paste0("Plots for the particlesInHeteroShape vs alpha and vs time"))

data <- read.table("chain_particlesInHeteroShape.tab", header=F)
data <- data[c(1,2,3)]
colnames(data) <- c("alpha","timestep","particlesInHeteroShape")
data$particlesInHeteroShape <- data$particlesInHeteroShape / totBeads
print(head(data))

colors <- rainbow(length(unique(data$alpha)))
#colors <- scale_colour_brewer(palette="Reds")(length(unique(data$alpha)))
print(paste0("Number of alphas ",length(unique(data$alpha))))

print(paste0("alphaShape"," ","volume"))
tolerance=0.005
outFileD <- "particlesInHeteroShape_vs_time.tab"
string <- data.frame("timestep","alphaPlateau","particlesInHeteroShapePlateau")
write.table(string,file=outFileD,row.names = F,quote = FALSE, col.names=F, sep="\t", append=T)

for(time in seq(0,100,1))
{
    # Plot the ashapeVol per replica for some values of timestep
    outFile <- paste0("chain_particlesInHeteroShape_vs_alpha_time_",time,".pdf")
    if(!file.exists(outFile))
    {
        print(paste0("outFile ",outFile))
	dataPlot <- data[data$timestep == time & data$alpha <= 15,]
	print(head(dataPlot))


	p <- ggplot(dataPlot, aes(x=factor(alpha), y=particlesInHeteroShape, fill="orange")) +
     	     geom_violin( position = position_dodge(1), trim=T, scale = "width") +
	     #scale_fill_manual(values=colors) +	     
	     geom_boxplot(position = position_dodge(1), fill="white", color="black", width=0.1, outlier.shape = NA) +     
	     ylim(0.0,1.0) +
	     theme(panel.background = element_rect(fill = NA),
                 panel.grid.major.x = element_blank(),
                 panel.grid.major.y = element_blank(),
            	 axis.title = element_text(face="bold",size=24),
		 axis.text = element_text(face="bold",size=12),
		 axis.text.x = element_text(angle=60, hjust=1),
            legend.position = "none") +
  	    labs(x="alpha",y="Fract. of particles in hetero shape",fill="") +     
	    guides(fill="none") #+
	    #geom_signif(comparisons = my_comparisons, y_position=c(ymax*(1+0.05),ymax*(1+0.05))) +
	    #geom_text(data=summ, aes(y=score, label = as.character(n)), color="black")

    	#p <- ggplot() +
	#    geom_violin(dataPlot, mapping=aes(x=factor(alpha), y=particlesInHeteroShape, color=rainbow(3)[[1]], alpha=0.2)) +	 
	#    scale_colour_manual(values =c(a=colorRampPalette(c('red',  'darkred'  ))(length(unique(data$replica))), b=colorRampPalette(c('green','darkgreen'))(length(un, c=colorRampPalette(c('blue',  'darkblue'))(length(unique(data$replica))))) +
        #    ylim(0.0,1.0) +
        #    theme(panel.background = element_rect(fill = NA),
        #         panel.grid.major.x = element_blank(),
        #         panel.grid.major.y = element_blank(),
        #    	 axis.title = element_text(face="bold",size=24),
	#	 axis.text = element_text(face="bold",size=18),
        #    legend.position = "none") +
        #    labs(x="alpha",y="Fract. particles in hetero shape",fill="") +     
        #    geom_text(data=stats,mapping=aes(x=alpha,y=0, label = n), color="black") +
        #    guides(fill="none")

    	pdf(outFile)
    	print(p)
    	dev.off()
    }
}

for(alpha in alphas)
{
    # Plot the ashapeVol per replica for some values of timestep
    outFile <- paste0("chain_particlesInHeteroShape_vs_time_alpha_",alpha,".pdf")
    if(!file.exists(outFile))
    {
        print(paste0("outFile ",outFile))
	dataPlot <- data[-1.0 < (data$alpha - alpha) & (data$alpha - alpha) < 0.1,]
	#dataPlot <- data[data$alpha == alpha,]
	print(head(dataPlot))

	p <- ggplot(dataPlot, aes(x=factor(timestep), y=particlesInHeteroShape, fill="orange")) +
     	     geom_violin( position = position_dodge(1), trim=T, scale = "width") +
	     #scale_fill_manual(values=colors) +	     
	     geom_boxplot(position = position_dodge(1), fill="white", color="black", width=0.1, outlier.shape = NA) +     
	     ylim(0.0,1.0) +
	     theme(panel.background = element_rect(fill = NA),
                 panel.grid.major.x = element_blank(),
                 panel.grid.major.y = element_blank(),
            	 axis.title = element_text(face="bold",size=24),
		 axis.text = element_text(face="bold",size=12),
		 axis.text.x = element_text(angle=60, hjust=1),
            legend.position = "none") +
  	    labs(x="time",y="Fract. of particles in hetero shape",fill="") +     
	    guides(fill="none") #+
	    #geom_signif(comparisons = my_comparisons, y_position=c(ymax*(1+0.05),ymax*(1+0.05))) +
	    #geom_text(data=summ, aes(y=score, label = as.character(n)), color="black")

    	#p <- ggplot() +
	#    geom_violin(dataPlot, mapping=aes(x=factor(timestep), y=particlesInHeteroShape, color=rainbow(3)[[1]], alpha=0.2)) +	 
	#    ylim(0.0,1.0) +
        #    theme(panel.background = element_rect(fill = NA),
        #         panel.grid.major.x = element_blank(),
        #         panel.grid.major.y = element_blank(),
        #    	 axis.title = element_text(face="bold",size=24),
	#	 axis.text = element_text(face="bold",size=18),
        #    legend.position = "none") +
        #    labs(x="time",y="Fract. particles in hetero shape",fill="") +     
            #geom_text(data=stats,mapping=aes(x=alpha,y=0, label = n), color="black") +
        #    guides(fill="none")

    	pdf(outFile)
    	print(p)
    	dev.off()
    }
}
quit()

# Compute stats on the dataset at a fixed timestep
stats <- data %>%
    group_by(timestep,alpha) %>%
    dplyr::reframe(n=n(), mean=mean(alphaShapeVol), sd=sd(alphaShapeVol), median=median(alphaShapeVol), mad=mad(alphaShapeVol))
stats <- unique(stats)
print(head(stats))

# Compute the derivative
t = 5000000
tmpData <- stats[stats$timestep == t,]
print(head(tmpData))
deriv <- data.frame(tmpData[1:length(tmpData$alpha)-1,]$alpha,diff(tmpData$mean)/diff(tmpData$alpha))
colnames(deriv) <- c("alpha","dMeanAlphaShapeVol")
print(head(deriv))
outFile <- paste0("chain_avg_alphaShape_dvol_vs_alpha_timestep_",t,".pdf")
if(!file.exists(outFile))
{
    print(paste0("outFile ",outFile))
    print(unique(data$replica))

    p <- ggplot() +
	 geom_line(deriv, mapping=aes(x=alpha, y=dMeanAlphaShapeVol)) + #, group=as.factor(timestep)), alpha=0.2) +
	 geom_line(tmpData, mapping=aes(x=alpha, y=mean)) + #, group=as.factor(timestep)), alpha=0.2) +	 
         ylim(0.0,1.0) +
	 #scale_colour_manual(values=rainbow(length(unique(data$timestep)))) +
         theme(panel.background = element_rect(fill = NA),
                 panel.grid.major.x = element_blank(),
                 panel.grid.major.y = element_blank(),
            	 axis.title = element_text(face="bold",size=24),
		 axis.text = element_text(face="bold",size=18),
         legend.position = "none") +
         labs(x="alpha",y="d(alpha-shape-volume)",fill="") +     
         guides(fill="none")

    pdf(outFile)
    print(p)
    dev.off()
}

quit()

# Compute the alpha end to zero point of the elbow point


# Compute the alpha end of the initial top plateau point of the elbow point


# Plot the ashapeVol averaged over replicates for all the values of timestep
outFile <- paste0("chain_avg_alphaShape_vol_vs_alpha.pdf")
if(!file.exists(outFile))
{
    print(paste0("outFile ",outFile))
    print(unique(data$replica))

    p <- ggplot() +
	 geom_line(data, mapping=aes(x=alpha, y=alphaShapeVol, group=as.factor(timestep)), alpha=0.2) +	 
         ylim(0.0,1.0) +
	 scale_colour_manual(values=rainbow(length(unique(data$timestep)))) +
         theme(panel.background = element_rect(fill = NA),
                 panel.grid.major.x = element_blank(),
                 panel.grid.major.y = element_blank(),
            	 axis.title = element_text(face="bold",size=24),
		 axis.text = element_text(face="bold",size=18),
         legend.position = "none") +
         labs(x="alpha",y="Fract. alpha-shape-volume",fill="") +     
         #geom_text(data=stats,mapping=aes(x=alpha,y=0, label = n), color="black") +
         guides(fill="none")

    pdf(outFile)
    print(p)
    dev.off()
}


quit()

#for(t in seq(0,5000000,5000))
#for(t in seq(0,5000000,250000))
for(t in unique(data$timestep))
{
    print(paste0("Time ",t))

    tmpData <- data[data$timestep == t,]
    print(head(tmpData))

    outFile <- paste0("chain_alphaShape_vol_vs_alpha_timestep_",t,"_perReplicate.pdf")
    if(!file.exists(outFile))
    {
        print(paste0("outFile ",outFile))
	print(head(tmpData))
	print(unique(tmpData$replica))

        p <- ggplot() +
     	     #geom_line(tmpData[tmpData$replica == "1A",], mapping=aes(x=alpha, y=alphaShapeVol, group=as.factor(replica), color=)) +
     	     geom_line(tmpData, mapping=aes(x=alpha, y=alphaShapeVol, group=as.factor(replica), color=replica), fill=rainbow(length(unique(data$replica)))) +	     
	     ylim(0.0,1.0) +
	     theme(panel.background = element_rect(fill = NA),
                 panel.grid.major.x = element_blank(),
                 panel.grid.major.y = element_blank(),
            	 axis.title = element_text(face="bold",size=24),
		 axis.text = element_text(face="bold",size=18),
           legend.position = "none") +
	   labs(x="alpha",y="Fract. alpha-shape-volume",fill="") +     
	   #geom_text(data=stats,mapping=aes(x=alpha,y=0, label = n), color="black") +
	   guides(fill="none")
	   
	   pdf(outFile) #, width=length(unique(tmpData$alpha))*2, height=14)
	   print(p)
	   dev.off()
    }
    next



    # Get the quantity and the alpha of the plateau    
    #print(tail(unique(stats),4))
    #quantity <- mean(tail(unique(stats),4)$mean)
    #print(quantity)

    #alphaT <- min(stats[abs(stats$mean-quantity)<tolerance,]$alpha)
    #string <- data.frame(t,alphaT,quantity)
    #write.table(string,file=outFileD,row.names = F,quote = FALSE, col.names=F, sep="\t", append=T)
    #print(paste0(t," ",alphaT," ",quantity))

    outFile <- paste0("chain_alphaShape_vol_vs_alpha_timestep_",t,".pdf")
    if(!file.exists(outFile))
    {
        print(paste0("outFile ",outFile))

        p <- ggplot(tmpData, aes(x=factor(alpha), y=alphaShapeVol, fill=factor(alpha), alpha=0.2)) +
             #geom_violin( tmpData, mapping=aes(x=factor(alpha), y=alphaShapeVol, fill=factor(alpha), alpha=0.2), position = position_dodge(1), trim=T, scale = "width") +
     	     geom_violin( position = position_dodge(1), trim=T, scale = "width") +     
     	     #geom_boxplot(tmpData, aes(x=factor(alpha), y=alphaShapeVol, fill=factor(alpha), color="black"), width=0.1, position = position_dodge(1), outlier.shape = NA) +
	     geom_boxplot(position = position_dodge(1), fill="white", color="black", width=0.1, outlier.shape = NA) +     
	     #scale_color_manual(fill="white") +
	     ylim(0.0,1.0) +
	     #scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
	     #labels = trans_format("log10", math_format(10^.x))) +
	     theme(panel.background = element_rect(fill = NA),
                 panel.grid.major.x = element_blank(),
                 panel.grid.major.y = element_blank(),
            	 axis.title = element_text(face="bold",size=24),
		 #axis.title = element_text(size=12),	   
		 axis.text = element_text(face="bold",size=18),
		 #axis.text.y = element_text(face="bold",size=18),
		 #axis.text.x = element_text(angle=60, hjust=1),
		 #axis.ticks.x = element_blank(),
		 #axis.ticks.y = element_blank(),
           legend.position = "none") +
	   labs(x="alpha",y="alpha-shape-volume / confining-sphere-volume",fill="") +     
	   #scale_color_manual(values=colors) +
	   #scale_fill_manual(values=colors) +
	   geom_text(data=stats,mapping=aes(x=as.factor(alpha),y=0, label = n), color="black") +
	   geom_line(data=stats,mapping=aes(x=as.factor(alpha),y=mean), color="black", group=1) +
	   geom_ribbon(data=stats,mapping=aes(x=as.factor(alpha),y = mean, ymin = mean - sd, ymax = mean + sd, fill="cyan", group=1), alpha = .5) +
	   geom_line(data=data.frame(x=as.factor(alpha),y=quantity),mapping=aes(x=x,y=y),color="black", linetype = "dashed")
	   geom_line(data=data.frame(x=as.factor(alphaT),y=c(0.0,1.0)),  mapping=aes(x=x,y=y),color="black", linetype = "dashed")	   
	   guides(fill="none")
	   
	   pdf(outFile, width=length(unique(tmpData$alpha))*2, height=14)
	   print(p)
	   dev.off()
     }

     quit()

} # Close cycle over t


print(paste0("Number of beads!"))
data <- read.table("chain_particlesInHeteroShape.tab", header=F)
data <- data[c(1,2,3,4)]
colnames(data) <- c("replica","timestep","alpha","alphaShapeBeads")
data$alphaShapeBeads <- data$alphaShapeBeads / totBeads
print(head(data))

tolerance=0.005
outFileD <- "alphaShapeBeads_vs_time.tab"
string <- data.frame("timestep","alphaPlateau","alphaShapeBeadsPlateau")
write.table(string,file=outFileD,row.names = F,quote = FALSE, col.names=F, sep="\t", append=T)

for(t in seq(0,5000000,250000))
#for(t in seq(0,5000000,5000))
#for(t in seq(0,5000000,5000000))
{
    #print(paste0("Time ",t))

    tmpData <- data[data$timestep == t,]
    #print(head(tmpData))

    # Compute stats on the dataset at a fixed timestep
    stats <- tmpData %>%
        group_by(alpha) %>%
        dplyr::reframe(x=factor(alpha), n=n(), mean=mean(alphaShapeBeads), sd=sd(alphaShapeBeads), median=median(alphaShapeBeads), mad=mad(alphaShapeBeads))
    #print(stats)

    # Get the quantity and the alpha of the plateau    
    #print(tail(unique(stats),4))
    quantity <- mean(tail(unique(stats),4)$mean)
    #print(quantity)

    alphaT <- min(stats[abs(stats$mean-quantity)<tolerance,]$alpha)
    string <- data.frame(t,alphaT,quantity)
    write.table(string,file=outFileD,row.names = F,quote = FALSE, col.names=F, sep="\t", append=T)

    #print(paste0(t," ",alphaT," ",quantity))
    next

    outFile <- paste0("chain_alphaShape_vol_vs_alpha_timestep_",t,".pdf")
    if(!file.exists(outFile))
    {
        print(paste0("outFile ",outFile))

        p <- ggplot(tmpData, aes(x=factor(alpha), y=alphaShapeBeads, fill=factor(alpha), alpha=0.2)) +
             #geom_violin( tmpData, mapping=aes(x=factor(alpha), y=alphaShapeBeads, fill=factor(alpha), alpha=0.2), position = position_dodge(1), trim=T, scale = "width") +
     	     geom_violin( position = position_dodge(1), trim=T, scale = "width") +     
     	     #geom_boxplot(tmpData, aes(x=factor(alpha), y=alphaShapeBeads, fill=factor(alpha), color="black"), width=0.1, position = position_dodge(1), outlier.shape = NA) +
	     geom_boxplot(position = position_dodge(1), fill="white", color="black", width=0.1, outlier.shape = NA) +     
	     #scale_color_manual(fill="white") +
	     ylim(0.0,1.0) +
	     #scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
	     #labels = trans_format("log10", math_format(10^.x))) +
	     theme(panel.background = element_rect(fill = NA),
                 panel.grid.major.x = element_blank(),
                 panel.grid.major.y = element_blank(),
            	 axis.title = element_text(face="bold",size=24),
		 #axis.title = element_text(size=12),	   
		 axis.text = element_text(face="bold",size=18),
		 #axis.text.y = element_text(face="bold",size=18),
		 #axis.text.x = element_text(angle=60, hjust=1),
		 #axis.ticks.x = element_blank(),
		 #axis.ticks.y = element_blank(),
           legend.position = "none") +
	   labs(x="alpha",y="alpha-shape--hetero-beads / total-beads",fill="") +     
	   #scale_color_manual(values=colors) +
	   #scale_fill_manual(values=colors) +
	   geom_text(data=stats,mapping=aes(x=as.factor(alpha),y=0, label = n), color="black") +
	   geom_line(data=stats,mapping=aes(x=as.factor(alpha),y=mean), color="black", group=1) +
	   geom_ribbon(data=stats,mapping=aes(x=as.factor(alpha),y = mean, ymin = mean - sd, ymax = mean + sd, fill="cyan", group=1), alpha = .5) +
	   #geom_line(data=data.frame(x=c(0,10.),y=c(quantity,quantity)),mapping=aes(x=x,y=y),color="black", linetype = "dashed")
	   #geom_line(data=data.frame(x=c(alphaT,alphaT),y=c(0.0,1.0)),  mapping=aes(x=x,y=y),color="black", linetype = "dashed")	   
	   guides(fill="none")
	   
	   pdf(outFile, width=length(unique(tmpData$alpha))*2, height=14)
	   print(p)
	   dev.off()
     }

     quit()

} # Close cycle over t
