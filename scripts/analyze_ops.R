library(psych)
library(xtable)
library(outliers)
library(tableplot)
library(ggplot2)
library(gridExtra)
library(reshape2)
#library(externalVector)
options(width=192)
CEXAXIS <- 0.51

args<-commandArgs(TRUE)
separator = "\n\n"
evoresults = read.table(args[1], header=TRUE)

#Prepare run data
fresults <- evoresults[ evoresults$Evaluations < 1000000 & !is.nan(evoresults$Best) & evoresults$Best < 0.01 & !is.na(evoresults$Problem),]
fresults$Problem <- factor(fresults$Problem)
fresults$Op <- factor(fresults$Op)
fresults$OpLength[fresults$Op == 'GCD'] <- 'other'
fresults$OpLength[fresults$Op == 'rnd7F'] <- 'other'
fresults$OpLength[fresults$Op != 'rnd7F' & fresults$Op != 'GCD'] <- fresults$OpSize[fresults$Op != 'rnd7F' & fresults$Op != 'GCD']
slabels = levels(fresults$OpSize)
slabels[1] = 'other'
slabels[6] = 'other'
slabels <- slabels[-1]
fresults$OpLength <- factor(fresults$OpLength,labels=slabels)
fresults$Rates <- factor(fresults$Rates)
fresults$kosher <- as.factor(paste(fresults$Op,fresults$Rates,fresults$OpLength,sep='.'))
temp = fresults[fresults$Op != 'rnd7F',]
splitres <- split(temp, temp$Problem)
#Prepare summary of data
sumresults = describeBy(fresults$Evaluations,interaction(fresults$Problem,fresults$Op,fresults$Rates, fresults$OpLength),mat=FALSE)
sumresultsb <- do.call("rbind",sumresults)
sumresultsb <- sumresultsb[order(rownames(sumresultsb)),]
rnames <- as.vector(rownames(sumresultsb))
splitted <- strsplit(rnames,c('[.]'))
splitm <- matrix(unlist(splitted), ncol=4, byrow=TRUE)
factorvars <- cbind(as.data.frame(splitm),sumresultsb)
#Filter problem with generalisation
genresults = rbind(splitres[["harmonic"]],splitres[[ "pendulum"]])
genresults = genresults[genresults$Validation < 1000000,]
genresults$Problem <- factor(genresults$Problem)
genresults$Op <- factor(genresults$Op)
genresults$OpLength <- factor(genresults$OpLength)
genresults$Rates <- factor(genresults$Rates)
genresults$kosher <- factor(genresults$kosher)

plotruns = '--plotruns' %in% args[]
if(plotruns){  
"PLOTRUNS
 Plots the boxplots of the number of evaluations and barplot with sucess rate for each problem."

splitres <- split(fresults, fresults$Problem)
for(i in 1:length(splitres)){
  prob = toString(splitres[[i]]$Problem[1])
  p <- ggplot(data = splitres[[i]][splitres[[i]]$Op != 'rnd7F',], aes(x = Rates, y = Evaluations,fill=OpLength)) + 
    geom_boxplot()+#outlier.shape=NA)+
    #coord_cartesian(ylim=c(0, 40000)) +
    facet_grid(~Op) +
    #scale_y_log10() +
    theme(axis.text.x=element_text(angle=-90))
  if(prob == 'pendulum')
    p <- p +  coord_cartesian(ylim=c(0, 100000))
  ggsave(paste(c('boxplot_', prob , '.pdf'),collapse=''))

  ggplot(data=splitres[[i]][splitres[[i]]$Op != 'rnd7F',], aes(x=Rates,fill=OpLength)) +
    geom_bar(stat='bin',position='dodge') +
     #coord_flip() +
      facet_grid(~Op) +
        theme(axis.text.x=element_text(angle=-90))
  ggsave(paste(c('success_', prob , '.pdf'),collapse=''))
}
}

success = '--success' %in% args[]
if(success){  
"PLOTSUCCESS
Plots the  barplot with success rate for each problem."

ggplot(data=fresults[fresults$Op != 'rnd7F',], aes(x=Rates,fill=OpLength)) +
  geom_bar(stat='bin',position='dodge') +
  #coord_flip() +
  facet_grid(Problem~Op) +
  theme(axis.text.x=element_text(angle=-90)) +
  ylab('Success Rate (%)')
ggsave('r1success.pdf')
}



compare = '--compare' %in% args[]
if(compare){  
"COMPARE
 Performs the statistical tests over the number of evaluations of the runs and plots the result in a graphical matrix."
print(factorvars)
textable = xtable(factorvars[c(1:4,6:9,12,17)],
                           digits=c(0, 0, 0, 0, 0, 0, 0, 0,0,0,0),
				  caption="Summary of the number of evaluations necessary to find an optimal solution for each experiment",
				  label="table:r0runs")
print(textable, include.rownames=FALSE, file = "summaryr1runs.tex")
cat(separator)


for(i in 1:length(splitres)){
  prob = toString(splitres[[i]]$Problem[1])
  print(prob)
  kr = kruskal.test(splitres[[i]]$Evaluations~splitres[[i]]$Op+splitres[[i]]$Rates+splitres[[i]]$OpLength)
  print(kr)
  #wb <- aggregate(splitres[[i]]$Evaluations,
             #   by = list(o = splitres[[i]]$Rates,
              #            r = splitres[[i]]$OpLength),
               # FUN = mean)
  #print(wb)
  #ft <- friedman.test(x ~ o | r, data = wb)
  #ft <- friedman.test(Evaluations ~ facts|Op, data = splitres[[i]])
  #print(ft)
  cat(separator)
  pw_kr = pairwise.wilcox.test(splitres[[i]]$Evaluations, splitres[[i]]$kosher,
    alternative='less',p.adj='bonferroni')
  pwtable <- melt(pw_kr[['p.value']])
  print(pw_kr)
  pw_kr2 = pairwise.wilcox.test(splitres[[i]]$Evaluations, splitres[[i]]$kosher,
    alternative='greater',p.adj='bonferroni')
  pwtable2 <- melt(pw_kr2[['p.value']])
  print(pw_kr2)
  #print(splitres[[i]]$kosher)
  n=length(unique(splitres[[i]]$kosher))
  pwtable$value2 <- pwtable2$value
  pwtable$Difference[pwtable$value <= 0.05] <- 'less'
  pwtable$Difference[pwtable$value2 <= 0.05] <- 'greater'
  pwtable$Difference[pwtable$value > 0.05 & pwtable$value2 > 0.05] <- 'indifferent'

  pwtable$Var1 <- factor(pwtable$Var1, levels=sort(unique(fresults$kosher), decreasing=TRUE))
  pwtable$Difference <- factor(pwtable$Difference, levels=c('greater','indifferent','less'))
  #pwtable$Var2 <- factor(pwtable$Var2, levels=sort(unique(factors), decreasing=FALSE))
  ggplot(data=pwtable, aes(x=Var2,y=Var1)) +
    geom_tile(aes(fill=Difference))+
    theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
    scale_fill_brewer(palette = 'PRGn',drop=FALSE) +
    theme(axis.text.x=element_text(angle=-90)) +
    xlab('') + ylab('') +
    opts(aspect.ratio = 1)
  ggsave(paste(c('r1mw_',prob,'.pdf'),collapse=''))
}
}

gen = '--gen' %in% args
if(gen){  
"GENERALISATION
 Performs the analysis of the Validation variable with statistical tests, and plots the results."
#print(genresults[[1]])
sumgen =  describeBy(genresults$Validation, interaction(genresults$Problem,genresults$Op,genresults$Rates, genresults$OpLength), mat=FALSE)

sumgenb <- do.call("rbind",sumgen)
sumgenb <- sumgenb[order(rownames(sumgenb)),]
rnames <- as.vector(rownames(sumgenb))
splitted <- strsplit(rnames,c('[.]'))
splitm <- matrix(unlist(splitted), ncol=4, byrow=TRUE)
factoredgen <- cbind(as.data.frame(splitm),sumgenb)
factoredgen$V1 <- factor(factoredgen$V1)
factoredgen$V2 <- factor(factoredgen$V2)
factoredgen$V3 <- factor(factoredgen$V3)
factoredgen$V4 <- factor(factoredgen$V4)

"
textable = xtable(sumgen[c(2,4:6,10)], 
                         digits=c(4, 4, 0, 5,5,6),
				  caption='Summary of the extrapolation results for each experiment.',
				  label='table:sumgen')
print(textable, file = 'summary-genkeijzer.tex')
"

spresults <- split(genresults,genresults$Problem)
#print(spresults)

for(i in 1:length(spresults)){
  #print(genresults[[i]]$Problem[1])
  #print(genresults[[i]])
  prob <- toString(spresults[[i]]$Problem[1])
  print(prob)
  ggplot(data = factoredgen[factoredgen$V1==prob,], aes(x=V3,y=mean,fill=V4)) +
    facet_grid(~V2) +
    #geom_line(aes(y=min)) +
      geom_point(aes(colour=V4),size=1.5,position=position_dodge(width=0.3)) +
        geom_errorbar(aes(ymin=min,ymax=max,colour=V4),width=.3,position='dodge') +
          theme(axis.text.x=element_text(angle=-90)) +
            #theme(legend.position = 'none') +
              xlab('Experiment') + ylab('Mean Fitness')
  ggsave(paste(c('r1gen_',prob,'.pdf'),collapse=''))

  if(prob == 'harmonic'){
    print(spresults[[i]][which(spresults[[i]]$Validation==min(spresults[[i]]$Validation)),])
  }
  else
    print(spresults[[i]][which(spresults[[i]]$Validation==max(spresults[[i]]$Validation)),])
  kr2 = kruskal.test(spresults[[i]]$Validation ~spresults[[i]]$kosher)
  print(kr2)
  cat(separator)
  #mann-whitney for each pair
  pw_kr = pairwise.wilcox.test(spresults[[i]]$Validation, spresults[[i]]$kosher, alternative='less', p.adj='bonferroni',paired=FALSE, na.rm=TRUE)
  print(pw_kr)
  pw_kr2 = pairwise.wilcox.test(spresults[[i]]$Validation, spresults[[i]]$kosher, alternative='greater', p.adj='bonferroni', paired=FALSE, na.rm=TRUE)
  print(pw_kr2)
  cat(separator)
  pwtable <- melt(pw_kr[['p.value']])
  pwtable2 <- melt(pw_kr2[['p.value']])
  pwtable$value2 <- pwtable2$value
  pwtable$Difference[pwtable$value < 0.05] <- 'less'
  pwtable$Difference[pwtable$value2 < 0.05] <- 'greater'
  pwtable$Difference[pwtable$value > 0.05 & pwtable$value2 > 0.05] <- 'indifferent'

  pwtable$Var1 <- factor(pwtable$Var1, levels=sort(unique(spresults[[i]]$kosher), decreasing=TRUE))
  pwtable$Difference <- factor(pwtable$Difference, levels=c('greater','indifferent','less'))
  #print(pwtable)
  #dev.new()
  ggplot(data=pwtable, aes(x=Var2,y=Var1)) +
    geom_tile(aes(fill=Difference))+
    theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
    scale_fill_brewer(palette = "PRGn" , drop = FALSE) +
    theme(axis.text.x=element_text(angle=-90)) +
    xlab('') + ylab('') +
    opts(aspect.ratio = 1)
  #prob <- toString(spresults[[i]]$Problem[1])
  ggsave(paste(c('r1mwgen_',prob,'.pdf'),collapse=''))
}
}

compare = '--baseline' %in% args[]
if(compare){  
"COMPARE BASELINE
 Performs the statistical tests over the number of evaluations of the runs versus the baseline results."
library(lawstat)
wt <- data.frame(problem = character(0), rname = character(0), pvalue = numeric(0))
for(i in 1:length(splitres)){
  prob = toString(splitres[[i]]$Problem[1])
  rnd7F <- fresults[fresults$Op == 'rnd7F' & fresults$Problem == prob,]
  #print(prob)
  allsplit <- split(splitres[[i]],splitres[[i]]$kosher,drop=TRUE)
  for(i in 1:length(allsplit)){
    #print(toString(unique(allsplit[[i]]$kosher)))
    wres <- wilcox.test(allsplit[[i]]$Evaluations,rnd7F$Evaluations,alternative='less')
    #data <- rbind(allsplit[[i]],rnd7F)
    #wres <- levene.test(data$Evaluations, data$Op,kruskal.test=TRUE)
    #print(wres)
    #wattr <- attributes(wres)
    #print(wattr)
    if(as.double(wres$p.value) <= 0.05/5){
      #if(as.integer(median(allsplit[[i]]$Evaluations)) < as.integer(median(rnd7F$Evaluations)))
        p.value = 'less'
      #else
       # p.value = 'greater'
    }
    else
      p.value = 'not significant'
    #p.value <- wres$p.value
    #print(prob)
    #print(unique(allsplit[[i]]$kosher))
    #print(p.value)
    #print(mean(allsplit[[i]]$Evaluations))
    #print(mean(rnd7F$Evaluations))
    wt = rbind(wt, as.data.frame(cbind(prob,p.value,toString(unique(allsplit[[i]]$kosher))))) 
  }
}
print(wt)
splitted <- strsplit(as.matrix(wt$V3),c('[.]'))
splitm <- matrix(unlist(splitted), ncol=3, byrow=TRUE)

wt2 <- data.frame(wt$prob,as.data.frame(splitm),as.matrix(wt$p.value),stringsAsFactors = FALSE)
colnames(wt2) <- c('V1','V2','V3','V4','pvalue')
#print(wt2)
#wt2$Result[wt2$pvalue < 0.05] <- 'less'#ifelse(wt$pvalue < 0.05,1,0)
#wt2$Result[wt2$pvalue >= 0.05] <- 'not significant'#
wt2$pvalue <- factor(wt2$pvalue, levels=c('greater','not significant','less'))
print(wt2)
ggplot(data=wt2, aes(x=V3,y=V4)) +
  facet_grid(V1~V2) +
    geom_tile(aes(fill=pvalue))+
    theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
    scale_fill_brewer(palette = 'PRGn',drop=FALSE) +
    theme(axis.text.x=element_text(angle=-90)) +
    xlab('') + ylab('') +
    opts(aspect.ratio = 1)
  #prob <- toString(spresults[[i]]$Problem[1])
  ggsave(paste(c('r1compbaseline.pdf'),collapse=''))
}

compare = '--test' %in% args
if(compare){  
"TESTING
 Performs the statistical tests over the number of evaluations of the runs versus the baseline results."
rnd7F <- fresults[fresults$Problem %in% genresults$Problem & fresults$Op == 'rnd7F',]
#rnd7F <- rnd7F[rnd7F$Op == 'rnd7F',]
print(rnd7F)
}

warnings()
quit()
