library(psych)
library(xtable)
library(outliers)
library(tableplot)
library(ggplot2)
library(doBy)
library(reshape2)
library(lawstat)
options(width=256)
CEXAXIS <- 0.51

args<-commandArgs(TRUE)
separator = "\n\n"
#Prepare runs data
evoresults = read.table(args[1], header=TRUE)
#evoresults = read.table("~/Documents/thesis/phdsupport/data/representation/analysis/grouped_evolution.txt", header=TRUE)
fresults = evoresults[ evoresults$Evaluations < 10000000 & evoresults$Problem != "<NA>" & !is.nan(evoresults$Best) & evoresults$Best < 0.01,]
#fresults$Problem[fresults$Problem == 'fibonacc-b',] = 'fibonacci'
#fresults$Problem[fresults$Problem == 'squares-rnd',] = 'squares-red'
fresults$Problem <- factor(fresults$Problem)
#fresults$Experiment <- factor(fresults$Experiment)
#Prepare summary stats
factors = list(fresults$Problem)
sumresults = describeBy(fresults$Evaluations, factors, na.rm=TRUE, mat=TRUE)
sumfuns = describeBy(fresults$NumFunctions, factors, mat=TRUE)#print(sumresults)
sumresults = merge(sumresults, sumfuns[c(2,5,10)], by='group1')
print(sumresults)
#sumproteins = describeBy(fresults$NumProteins, factors, mat=TRUE)#print(sumresults)
#sumresults = merge(sumresults, sumproteins[c(2,5)], by='group1')
#print(sumresults)
#ordersum <- orderBy(~group1, data=sumresults)
#Filter problems with generalisation
spresults <- split(fresults,fresults$Problem)
genresults <- fresults
#genresults <- genresults[-which((genresults$Problem == 'nbitparity' | genresults$Problem=='squares-rnd' | genresults$Problem=='squares-full' | genresults$Problem=='modfactorial-k2s1-red' | genresults$Problem=='modfactorial-k2s2-red' | genresults$Problem=='modfactorial-k2s3-red' )& genresults$Validation > 0),]
#genresults = genresults[genresults$Validation < 1000000,]
genresults <- genresults[genresults$Validation == 0,]
genresults$Problem <- factor(genresults$Problem)

#genresults$Experiment <- factor(genresults$Experiment)
#Prepare validation stats
sumgen =  describeBy(genresults$Validation, genresults$Problem, mat=TRUE)
sumgen$success <- sumgen$n / sumresults$n
print(sumgen)

plotruns = '--plotruns' %in% args[]
if(plotruns){  
"PLOTRUNS
 Plots the boxplots of the number of evaluations and barplot with sucess rate for each problem, and the error bars for the #evaluations."
ggplot(data = fresults, aes(x = Problem, y = Evaluations)) + 
  geom_boxplot()+
  #facet_grid(Problem ~ .) +
  theme(axis.text.x=element_text(angle=-90))+
  theme(legend.position = "none")
ggsave('boxplot_evaluations.pdf')

ggplot(data = sumresults, aes(x=group1, y = mean.x, ymin=mean.x-se,ymax=mean.x+se)) +
  #facet_grid(group1 ~ .) +
  geom_line() +
  geom_point() +
  geom_errorbar(width=.3,position='dodge') +
  theme(axis.text.x=element_text(angle=-90)) +
  theme(legend.position = "none") +
  xlab('Experiment') +
  ylab('Mean number of evaluations')
ggsave('errorbars_evaluations.pdf')

print(sumresults)
textable = xtable(sumresults[c(1,4:6,10,15,16,17)],
                           digits=c(0, 0, 0, 0, 0, 0, 0, 1, 0),
				  caption="Summary of the number of evaluations necessary to find an optimal solution for each experiment",
				  label="table:r0runs")
print(textable, include.rownames=FALSE, file = "summaryruns.tex")
cat(separator)

ggplot(data = sumresults, aes(x = group1, y=n)) + 
  geom_bar(stat='identity')+
  #facet_grid(group1 ~ .)+
  theme(axis.text.x=element_text(angle=-90)) +
  theme(legend.position = "none") +
  xlab('Experiment') +
  ylab('Success Rate (%)')
ggsave('success.pdf')
}


compare = '--compare' %in% args[]
if(compare){  
"COMPARE
 Performs the statistical tests over the number of evaluations of the runs and plots the result in a graphical matrix."

t <- wilcox.test(spresults[['squares-full']]$Evaluations,spresults[['squares-red']]$Evaluations, alternative='less')
print(t)

}

gen = '--gen' %in% args
if(gen){  
"GENERALISATION
 Performs the analysis of the Validation variable with statistical tests, and plots the results."

textable = xtable(sumgen[c(2,16)], 
                         digits=c(0,4, 2),
				  caption="Summary of the extrapolation results for each experiment in the harmonic and pendulum problems.",
				  label="table:sumgen")
print(textable, include.rownames=FALSE,file = "summarygen.tex")
#print(sumgen)

ggplot(data = sumgen, aes(x = group1, y=success)) + 
  geom_bar(stat='identity')+
  #facet_grid(group1 ~ .)+
  theme(axis.text.x=element_text(angle=-90)) +
  theme(legend.position = "none") +
  xlab('Experiment') +
  ylab('Generalisation Success Rate (%)')
ggsave('gensuccess.pdf')

}


compare = '--mins' %in% args
if(compare){
  r <- data.frame()
  sp <- split(genresults, genresults$Problem)

  for(i in 1:length(sp)){
    r <- rbind(r,sp[[i]][sp[[i]]$NumFunctions == min(sp[[i]]$NumFunctions),][1,])
    #print()
  }
  print(r)
}



compare = '--test' %in% args
if(compare){  
print("TESTING")
prob = 'pendulum'
#prob = 'harmonic'
#np <- fresults$NumProteins
#dev.new()

#l <- levene.test(fresults$Evaluations,fresults$Experiment,kruskal.test=TRUE)
#print(l)
#print(fresults[1:10])

#fresults$Evaluations <- factor(fresults$Evaluations)
#fresults <- fresults[order(fresults$Evaluations),]
#plot(fresults[fresults$Problem=='harmonic',]$Evaluations)
#wt <- data.frame(problem=character(0),exp2=character(0),exp1=character(0),pvalue=numeric(0),median2=numeric(0),median1=numeric(0))
"
wt <- data.frame()
for(i in 1:length(spresults)){
  exp <- split(spresults[[i]],spresults[[i]]$Experiment)
  prob <- toString(spresults[[i]]$Problem[1])
  for(j in 1:length(exp)){
    for(k in 1:length(exp)){
      if( j < k ){
        data <- rbind(exp[[j]],exp[[k]])
        l <- levene.test(data$Evaluations,data$Experiment,kruskal.test=TRUE)
        wt <- rbind(wt, cbind(problem = prob, exp2 = toString(exp[[j]]$Experiment[1]),exp1 = toString(exp[[k]]$Experiment[1]),pvalue = l[['p.value']],median2 = median(exp[[j]]$Evaluations),median1 = median(exp[[k]]$Evaluations)))
      }
    } 
  }
}

#print(wt)
wtsp <- split(wt, wt$problem)
ncomp <- as.numeric(length(wt$problem[wt$problem == 'harmonic']))
print(paste(c('Number of pairwise comp.: ', ncomp),collapse=''))
#wt$median1 <- as.numeric(wt$median1)
#wt$median2 <- as.numeric(wt$median2)
print(0.05/ncomp)
lala <- as.vector(0.05/ncomp)
wt$pvalue <- as.numeric(levels(wt$pvalue)[as.double(wt$pvalue)])
wt$median1 <- as.numeric(levels(wt$median1)[as.integer(wt$median1)])
wt$median2 <- as.numeric(levels(wt$median2)[as.integer(wt$median2)])
#print(as.numeric(levels(wt$pvalue)[as.double(wt$pvalue)]) < lala)
#print(as.vector(wt$median1) < as.vector(wt$median2))

wt$Difference[wt$pvalue < 0.05/ncomp & wt$median1 < wt$median2] <- 'less'
wt$Difference[wt$pvalue < 0.05/ncomp & wt$median1 > wt$median2] <- 'greater'
wt$Difference[wt$pvalue > 0.05/ncomp ] <- 'indifferent'
print(wt)
wtsp <- split(wt, wt$problem)
for(i in 1:length(wtsp)){
  prob <- toString(wtsp[[i]]$problem[1])
  #wtsp[[i]] <- wtsp[[i]][orderBy(wtsp[[i]]$exp1,reverse=TRUE),]
  ggplot(data=wtsp[[i]], aes(x=exp2,y=exp1)) +
    geom_tile(aes(fill=Difference))+
    theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
    scale_fill_brewer(palette = 'PRGn'',drop=FALSE) +
    theme(axis.text.x=element_text(angle=-90)) +
    xlab('') + ylab('') +
    opts(aspect.ratio = 1)
  #prob <- toString(spresults[[i]]$Problem[1])
  ggsave(paste(c('r0brown_',prob,'.pdf'),collapse=''))
}
"
for(i in 1:length(spresults)){
  prob <- toString(spresults[[i]]$Problem[1])
  print(spresults[[i]]$Problem[1])
  exp <- split(spresults[[i]],spresults[[i]]$Experiment)
  for(j in 1:length(exp)){
    print(exp[[j]]$Experiment[1])
    ggplot(data = exp[[j]], aes(x = Evaluations)) + 
      geom_histogram()+
          theme(axis.text.x=element_text(angle=-90)) +
            scale_x_log10()+
          theme(legend.position = 'none')# +
#  xlab('Experiment') +
#  ylab('Success Rate (%)')
    
    ggsave(paste(c('evals_', i, j, '_', prob, '.pdf'),collapse=''))
    if(nrow(exp[[j]])>2){
      sw <- shapiro.test(log10(exp[[j]]$Evaluations))
      print(sw)
    }
  }
}

"
ggplot(genresults[genresults$Problem == prob,], aes(x=Experiment,y=Validation)) +
  geom_point(aes(colour=NumFunctions)) +
    scale_color_gradient2(low='blue',high='red')
    #geom_smooth(method='lm') +
      #scale_x_log10() +
        #scale_y_log10()
ggsave('test.pdf')

harm <- sumresults[sumresults$group1 == prob,]
harmgen <- sumgen[sumgen$group1 == prob,]
harm$mean.validation <- hqarmgen$mean
#print(sumgen)
#print(harm)
dev.new()
ggplot(harm, aes(x=mean,y=mean.validation)) +
  geom_point() +
    #geom_smooth(method='lm') +
      scale_x_log10() +
        scale_y_log10()
}
"
}
warnings()
traceback()
quit()

"
print('###MUTATIONS###')
##
mutresults = fresults[fresults$NeutralMut > 0,]
#print(mutresults)
mutsum = describeBy(mutresults$NeutralMut, mutresults$Experiment, mat=TRUE, na.rm = TRUE)
bitsum = describeBy(mutresults$NeutralBits, mutresults$Experiment, mat=TRUE, na.rm = TRUE)
sizesum =  describeBy(mutresults$AvgGeneSize, mutresults$Experiment, mat=TRUE, na.rm = TRUE)
print(mutsum)
#dev.new()
#error.bars(stats=mutsum, labels=mutsum$group1, ylab='Neutral Mutations Rate' )

dev.new()
plot(mutsum$group1, mutsum$mean)
points(bitsum$group1, bitsum$mean)


print('###PROTEINS###')
sumpnum = describeBy(fresults$AvgNumProteins, fresults$Experiment, mat=TRUE)
sumfnum = describeBy(fresults$AvgNumFunctions, fresults$Experiment, mat=TRUE)
dev.new()
counts = data.frame(sumpnum$mean,sumfnum$mean)
print(counts)
barplot(t(counts), names.arg=sumpnum$group1,ylab='#Proteins/Functions', xlab='Experiment', main = 'Number of Proteins and Functions in the Genomes', las=3, col=c('darkblue','red'), legend = list('#Proteins','#Functions'), beside=TRUE, ylim=c(0,35))

dev.new()
plot(sumpnum$mean,sumresults$mean, ylab = '#Evaluations', xlab='#Proteins')

dev.new()
plot(sumpnum$group1, sumfnum$mean/sumpnum$mean,ylab='#Functions/#Proteins')

#dev.new()
#plot((sumfnum$mean/sumpnum$mean),sumgen$mean)

dev.new()
plot(sumpnum$mean, sumgen$mean, ylab='#Proteins', xlab='Generalization Error')
dev.new()
plot(sumpnum$mean, sumgen$n)
#plot(sumpnum$mean, sumresults$mean)
#error.bars(stats=sumpnum, label=sumpnum$group1, add=TRUE)
#dev.new()
#contour(mutsum$mean, bitsum$mean, sizesum$mean )
"
