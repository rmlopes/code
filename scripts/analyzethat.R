library(psych)
library(xtable)
library(outliers)
library(tableplot)
library(ggplot2)
library(doBy)
options(width=256)
CEXAXIS <- 0.51

args<-commandArgs(TRUE)
separator = "\n\n"
#evoresults = read.table(args[1], header=TRUE)
evoresults = read.table("~/Documents/thesis/phdsupport/data/representation/analysis/grouped_evolution.txt", header=TRUE)
fresults = evoresults[ evoresults$Evaluations < 1000000 & evoresults$Problem != "<NA>" & !is.nan(evoresults$Best) & evoresults$Best < 0.01,]
#|| evoresults$Problem == 'harmonic',]
#fresults = fresults[fresults$Best < 0.01,]
#print(fresults[0:10])
#fresults <- rbind(evoresults[evoresults$Problem == "harmonic",], fresults)
fresults$Problem <- factor(fresults$Problem)
fresults$Experiment <- factor(fresults$Experiment)
dev.new()
#b = boxplot(Evaluations ~interaction(Problem,Experiment), fresults,las=3, outline = FALSE, main="Number of Evaluations to succeed",col=(c("green","red",'blue','yellow')))
ggplot(data = fresults, aes(x = Experiment, y = Evaluations, fill = Problem)) + 
  geom_boxplot()+
  facet_grid(Problem ~ .) +
  theme(axis.text.x=element_text(angle=-90))+
  theme(legend.position = "none")

factors = list(fresults$Problem, fresults$Experiment)
sumresults = describeBy(fresults$Evaluations, factors, na.rm=TRUE, mat=TRUE)
print(sumresults)

dev.new()
ggplot(data = sumresults, aes(x=group2, y = mean, colour=group1, ymin=mean-se,ymax=mean+se)) +
  facet_grid(group1 ~ .) +
  geom_line() +
  geom_point() +
  geom_errorbar(width=.3,position='dodge') +
  theme(axis.text.x=element_text(angle=-90)) +
  theme(legend.position = "none")
#error.bars(fresults$Evaluations, stats=sumresults, xlab=sumresults$group1+":"+sumresults$group2, ylab="Number of Evaluations with 95% confidence" )
ordersum <- orderBy(~group1+group2, data=sumresults)
textable = xtable(ordersum[c(2,3,5:7,10)],
                           digits=c(0, 0, 0, 0, 0, 0, 0),
				  caption="Summary of the results for the Symbolic Regression experiments.",
				  label="table:summary-symb")
print(textable, file = "summary-symb.tex")
cat(separator)
normcount = sumresults$n /100
dev.new()
#barplot(normcount, main='Successful Runs(%)', names.arg = sumresults$group1, cex.names=CEXAXIS)
ggplot(data = sumresults, aes(x = group2, y=n/100, fill=group1)) + 
  geom_bar(stat='identity')+
  facet_grid(group1 ~ .)+
  theme(axis.text.x=element_text(angle=-90)) +
  theme(legend.position = "none")

#for univariate normality testing
#qqnorm(fresults$V7[V1 == 'symb'])
#qqline(fresults$V7[V1 == 'symb'])

#non parametric
spresults <- split(fresults,fresults$Problem)
#spresults$Experiment <- factor(spresults$Experiment)
for (i in 1:length(spresults)) {
  print(spresults[[i]]$Problem[1])
  kr = kruskal.test(spresults[[i]]$Evaluations ~spresults[[i]]$Experiment)
  print(kr)
  #mann-whitney for each pair
  pw_kr = pairwise.wilcox.test(spresults[[i]]$Evaluations, spresults[[i]]$Experiment,
  alternative="less",p.adj="bonferroni",paired=FALSE,exact=FALSE)
  print(pw_kr)
  dev.new()
  pmatrix = pw_kr[['p.value']]
  pmatrix <- pmatrix[order(rownames(pmatrix))]
  heatmap(pw_kr[['p.value']])
  dev.off()
  pw_kr2 = pairwise.wilcox.test(spresults[[i]]$Evaluations, spresults[[i]]$Experiment,
  alternative="greater",p.adj="bonferroni",paired=FALSE,exact=FALSE)
  print(pw_kr2)
  #zz = qnorm(pw_kr$p.value)
  #print(zz)
  print("parametric tests...")
  aovresults <- aov(Evaluations~Experiment,data=spresults[[i]])
  print(summary(aovresults))
  tuk = TukeyHSD(aovresults)
  print(tuk)
  #dev.new()
  #plot(tuk)
  spnames <- strsplit(rownames(tuk$Experiment),'-')
  spnames0 <- sapply(spnames,'[', 1)
  spnames1 <- sapply(spnames,'[', 2)
  newtuk <- data.frame(as.numeric(tuk$Experiment[,4]),spnames0,spnames1)
  #sptuk <- rbind(tuk$Experiment,spnames)
  #newtuk <- as.matrix(newtuk)
  colnames(newtuk) <- c('padj','x','y')
  newtuk$x <- factor(newtuk$x)
  newtuk$y <- factor(newtuk$y)
  #print(newtuk)
  newtuk$padj1 <- cut(newtuk$padj,breaks=c(0,0.05,1.1),right=FALSE)
 # dev.new()
  ggplot(data=newtuk, aes(x, y), colour='white') +
    geom_tile(aes(fill=padj1)) +
      scale_fill_brewer(palette = "PRGn")+
        theme(axis.text.x=element_text(angle=-90)) +
          labs(x='', y='')
  ggsave(paste(c('tukey',i,'.pdf'),collapse=''))
 
                                        # +
      #scale_fill_gradient(low = "white",high = "steelblue")
  #heatmap(as.matrix(newtuk), Rowv=NA, Colv=NA, col = cm.colors(256), scale="column", margins=c(5,10))
  #cat(separator)
  
}

#quit()
#
#parametric
#ow = oneway.test(fresults$Evaluations~fresults$Experiment)
#print(ow)#pairwise; what about correction?
#cat(separator)
#changing p.adj results in apparently incoherent p-values
#pw_ow = pairwise.t.test(fresults$Evaluations, fresults$Experiment, p.adj="none",pool.sd=F)
#print(pw_ow)
#cat(separator)

##GENERALIZATION
print('###GENERALIZATION###')
#genresults = rbind(fresults[fresults$Validation > 0.0,] , evoresults[evoresults$Problem == 'harmonic',])
#genresults = rbind(fresults[fresults$Problem == "harmonic",], fresults[fresults$Problem == "pendulum",])
genresults = rbind(spresults[["harmonic"]],spresults[[ "pendulum"]])
genresults = genresults[genresults$Validation < 1000000,]
#genresults = fresults
genresults$Problem <- factor(genresults$Problem)
genresults$Experiment <- factor(genresults$Experiment)
sumgen =  describeBy(genresults$Validation, list(genresults$Problem,genresults$Experiment), mat=TRUE)
textable = xtable(sumgen[c(2,4:6,10)], 
                         digits=c(4, 4, 0, 5,5,6),
				  caption="Summary of the extrapolation results for each experiment.",
				  label="table:sumgen")
print(textable, file = "summary-genkeijzer.tex")
#print(sumgen)
pgen <- sumgen[sumgen$group1 == "pendulum",]
pgen$mean <- pgen$mean / max(pgen$mean)
pgen$se <- pgen$se / max(pgen$se)
hgen <- sumgen[sumgen$group1 == "harmonic",]
hgen$mean <- hgen$mean / max(hgen$mean)
hgen$se <- hgen$se / max(hgen$se,na.rm=TRUE)
normgen <- rbind(pgen,hgen)
#sumgen <- normalise(sumgen, "mean", by="group1")
#print(normgen)
#normcount = sumgen$n / sumresults$n
#dev.new()
#barplot(normcount, main='Success(%)', newtuk.arg = sumresults$group1, cex.newtuk=CEXAXIS)

#dev.new()
#error.bars(stats=sumgen, labels=sumgen$group2, ylab="Generalization Error" )
#print(normgen)
dev.new()
ggplot(data = normgen, aes(x=group2, y = mean, colour=group1, ymin=mean-se,ymax=mean+se)) +
  facet_grid(group1 ~ .) +
  geom_line() +
  geom_point() +
  geom_errorbar(width=.3,position='dodge') +
  theme(axis.text.x=element_text(angle=-90)) +
  theme(legend.position = "none")

genresults <- split(genresults,genresults$Problem)
for(i in 1:length(genresults)){
  print(genresults[[i]]$Problem[1])
  kr2 = kruskal.test(genresults[[i]]$Validation ~genresults[[i]]$Experiment)
  print(kr2)
  cat(separator)
  #mann-whitney for each pair
  pw_kr2 = pairwise.wilcox.test(genresults[[i]]$Validation, genresults[[i]]$Experiment, alternative="less", p.adj="bonferroni",paired=FALSE)
  print(pw_kr2)
  pw_kr2 = pairwise.wilcox.test(genresults[[i]]$Validation, genresults[[i]]$Experiment, alternative="greater", p.adj="bonferroni", paired=FALSE)
  print(pw_kr2)
  cat(separator)
}


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
