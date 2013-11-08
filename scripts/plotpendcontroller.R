library(psych)
library(xtable)
library(outliers)
library(tableplot)
library(ggplot2)
library(doBy)
library(reshape2)
library(grid)
library(gridExtra)
library(gtable)
options(width=256)
CEXAXIS <- 0.51

#grid.newpage()
args<-commandArgs(TRUE)
separator = "\n\n"
#evoresults = read.table(args[1], header=TRUE)
validation = as.matrix(read.table(args[1], header=FALSE))
vars = read.table(args[2], header=FALSE)
tval = t(validation)
all = cbind(vars,stack(as.data.frame(tval)))
all$V1 <- factor(all$V1)
all$V2 <- factor(all$V2,levels=sort(unique(all$V2),decreasing=TRUE))
all$V3 <- factor(all$V3)#,levels=sort(unique(all$V3),decreasing=TRUE))
all$V4 <- factor(all$V4,levels=sort(unique(all$V4),decreasing=FALSE))
all$values <- factor(all$values)
all$Result <- ifelse(all$values==0,"Sucessful","Failed")
print(all)
mylabels = c('Angular velocity (º/s)','Velocity (m/s)',"Cart Position (m)","Pole Angle (rad)")

fourvarplot <- function(alldata,labels){
  "Input must have the variables named V1,V2,V3,V4, and categorical Result"
pleg <- ggplot(data=alldata, aes(x=V4,y=V2)) +
  facet_grid(V1 ~ V3) + 
  geom_tile(aes(fill=Result))+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
  scale_fill_brewer(palette = "PRGn") +
  theme(axis.text.x=element_text(angle=-90)) +
  xlab(labels[1]) + ylab(labels[2]) +
  opts(aspect.ratio = 1) +
  theme(legend.position='right')

leg = gtable_filter(ggplot_gtable(ggplot_build(pleg)), "guide-box")

p <-  ggplot(data=alldata, aes(x=V4,y=V2)) +
  facet_grid(V1 ~ V3) + 
  geom_tile(aes(fill=Result))+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
  scale_fill_brewer(palette = "PRGn") +
  theme(axis.text.x=element_text(angle=-90)) +
  xlab(labels[1]) + ylab(labels[2]) +
  opts(aspect.ratio = 1) +
  theme(legend.position='none')

# get gtable object
z <- ggplot_gtable(ggplot_build(p))

# add label for right strip
z <- gtable_add_cols(z, z$widths[[8]])
z <- gtable_add_grob(z, 
  list(rectGrob(gp = gpar(col = NA, fill = gray(0.8))),
  textGrob(labels[3], rot = -90, gp = gpar(col = gray(1)))),
  4, 14, 12, 14, name = paste(runif(2)))

# add label for top strip
z <- gtable_add_rows(z, z$heights[[3]], 2)
z <- gtable_add_grob(z, 
  list(rectGrob(gp = gpar(col = NA, fill = gray(0.8))),
  textGrob(labels[4], gp = gpar(col = gray(1)))),
  3, 4, 3, 12, name = paste(runif(2)))

#add the legend now
z <- gtable_add_cols(z, z$widths[[8]])
z <- gtable_add_grob(z,
                     leg,
                     4, 15, 12, 16)

# add margins
z <- gtable_add_cols(z, unit(1/8, "line"), 8)
z <- gtable_add_rows(z, unit(1/8, "line"), 3)
return(z)
}

grid.newpage()
grid.draw(fourvarplot(all,mylabels))
#
#ggsave('r0bestpendgen.pdf',plot=p)
quit()
