options(scipen=9999)
library(dplyr)
library(ggplot2)
library(ggpubr)
library("ggrepel")
library(scales)
library("RColorBrewer")
args = commandArgs(trailingOnly=TRUE)
pi=3.14159265359

ymax <- as.numeric(args[[1]])

data <- read.table("cisCounts.tab", header=T)
print(head(data))
#quit()
print(unique(data$time))

levels_sim <- unique(data$time)

colors <- c(rep(c("white"),as.integer(length(levels_sim))))

levels <- levels_sim
print(levels)
print(length(levels))
print(length(colors))

data$time<- factor(data$time, levels = levels)
data <- data[!is.na(data$time),]

mC <- levels
for(time in mC)
{
    #print(time)
    result <- c(time,as.numeric(quantile(data[data$time==time,]$cisCounts, c(.25, .50, .75), na.rm = TRUE)))
    #print(result)

    if (!exists("quartiles"))
    {
	print(paste0("Adding ",time))
	quartiles <- data.frame()    
        quartiles <- as.data.frame(t(result))
    } else {
	print(paste0("Adding ",time))
        quartiles <- rbind(quartiles,as.data.frame(t(result)))
    }
}
colnames(quartiles) <- c("time","first","second","third")
print(quartiles)

for(time in unique(data$time))
{
    print(time)
    result <- c(time,as.numeric(quantile(data[data$time==time,]$cisCounts, c(.25, .50, .75), na.rm = TRUE)))
    print(result)
}
#quit()

# Plot the ashapeVol per replica for some values of timestep
outFile <- paste0("cisCounts.pdf")
if(!file.exists(outFile))
{
    print(paste0("outFile ",outFile))

    p <- ggplot(data, aes(x=factor(time), y=cisCounts, fill=factor(time))) +
             #DMSOannotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = as.numeric(quartiles[c(1),]$first), ymax = as.numeric(quartiles[c(1),]$third), fill = rgb(235/255,238/255,240/255,0.80), color=NA) + # DMSO
	     #DMSOannotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = as.numeric(quartiles[c(2),]$first), ymax = as.numeric(quartiles[c(2),]$third), fill = rgb(1,0,0,0.30), color=NA) + # DMSO_A
	     #DMSOannotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = as.numeric(quartiles[c(3),]$first), ymax = as.numeric(quartiles[c(3),]$third), fill = rgb(0,0,1,0.30), color=NA) + # DMSO_B
             #TSAannotate(geom = "rect", ymin = as.numeric(quartiles[c(4),]$first), ymax = as.numeric(quartiles[c(4),]$third), xmin = -Inf, xmax = Inf, fill = rgb(233/255,200/255,173/255,0.80), color=NA) + # TSA
             #TSAannotate(geom = "rect", ymin = as.numeric(quartiles[c(5),]$first), ymax = as.numeric(quartiles[c(5),]$third), xmin = -Inf, xmax = Inf, fill = rgb(1,0,0,0.30), color=NA) + # TSA_A
             #TSAannotate(geom = "rect", ymin = as.numeric(quartiles[c(6),]$first), ymax = as.numeric(quartiles[c(6),]$third), xmin = -Inf, xmax = Inf, fill = rgb(0,0,1,0.30), color=NA) + # TSA_B
     	     geom_violin( position = position_dodge(1), trim=T, scale = "width") +
	     scale_fill_manual(values=colors) +	     
	     geom_boxplot(position = position_dodge(1), fill="white", color="black", width=0.1, outlier.shape = NA) +     
	     ylim(0,ymax) +
	     #scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
	     #labels = trans_format("log10", math_format(10^.x))) +
	     theme(panel.background = element_rect(fill = NA),
                 panel.grid.major.x = element_blank(),
                 panel.grid.major.y = element_blank(),
            	 axis.title = element_text(face="bold",size=24),
		 #axis.title = element_text(size=12),	   
		 axis.text = element_text(face="bold",size=12),
		 #axis.text.y = element_text(face="bold",size=18),
		 axis.text.x = element_text(angle=60, hjust=1),
		 #axis.ticks.x = element_blank(),
		 #axis.ticks.y = element_blank(),
           legend.position = "none") +
	   labs(x="",y="cis-chromosome counts",fill="") +     
	   guides(fill="none")

    #pdf(outFile,width=28)
    pdf(outFile)    
    print(p)
    dev.off()

}

quit()
