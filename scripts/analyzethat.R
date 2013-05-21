library(psych)
library(xtable)
library(outliers)
library(tableplot)
options(width=192)
CEXAXIS <- 0.51

args<-commandArgs(TRUE)
separator = "\n\n"
evoresults = read.table(args[1], header=TRUE)
#evoresults = read.table("~/Documents/thesis/thedvd/experimentaldata/rencode/representation/keijzer_evolution.txt", header=TRUE)
fresults = evoresults[evoresults$Evaluations < 1000000,]
fresults = fresults[fresults$Best < 0.001,]
print(fresults[0,])
dev.new()
b = boxplot(Evaluations ~Experiment, fresults,cex.axis=CEXAXIS, outline = TRUE, main="Number of Evaluations to succeed")

sumresults = describeBy(fresults$Evaluations, fresults$Experiment, mat=TRUE)
dev.new()
error.bars(fresults$Evaluations, stats=sumresults, labels=sumresults$group1, ylab="Evaluations" )

textable = xtable(sumresults[c(2,4:6,10)],
                           digits=c(0, 0, 0, 0,0,0),
				  caption="Summary of the results for the Symbolic Regression experiments.",
				  label="table:summary-symb")
print(textable, file = "summary-symb.tex")
print(sumresults[c(2,4:6,10)])
cat(separator)
normcount = sumresults$n /100
dev.new()
barplot(normcount, main='Successful Runs(%)', names.arg = sumresults$group1, cex.names=CEXAXIS)

#for univariate normality testing
#qqnorm(fresults$V7[V1 == 'symb'])
#qqline(fresults$V7[V1 == 'symb'])

#non parametric
kr = kruskal.test(fresults$Evaluations ~fresults$Experiment)
print(kr)
cat(separator)
#mann-whitney for each pair
pw_kr = pairwise.wilcox.test(fresults$Evaluations, fresults$Experiment,
  alternative="less",p.adj="bonferroni")
print(pw_kr)
pw_kr2 = pairwise.wilcox.test(fresults$Evaluations, fresults$Experiment,
  alternative="greater",p.adj="bonferroni")
print(pw_kr2)
#zz = qnorm(pw_kr$p.value)
#print(zz)
cat(separator)

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
genresults = fresults[fresults$Generalization < 0.001,]
#genresults = fresults
sumgen =  describeBy(genresults$Generalization, genresults$Experiment, mat=TRUE)
textable = xtable(sumgen[c(2,4:6,10)], 
                         digits=c(4, 4, 0, 5,5,6),
				  caption="Summary of the extrapolation results for each experiment.",
				  label="table:sumgen")
print(textable, file = "summary-genkeijzer.tex")
print(sumgen)

#normcount = sumgen$n / sumresults$n
#dev.new()
#barplot(normcount, main='Success(%)', names.arg = sumresults$group1, cex.names=CEXAXIS)

dev.new()
error.bars(stats=sumgen, labels=sumgen$group1, ylab="Generalization Error" )


kr2 = kruskal.test(genresults$Generalization ~genresults$Experiment)
print(kr2)
cat(separator)
#mann-whitney for each pair
pw_kr2 = pairwise.wilcox.test(genresults$Generalization, genresults$Experiment, alternative="less", p.adj="bonferroni")
print(pw_kr2)
cat(separator)



print('###MUTATIONS###')
##
mutresults = fresults[fresults$NeutralMut > 0,]
#print(mutresults)
mutsum = describeBy(mutresults$NeutralMut, mutresults$Experiment, mat=TRUE, na.rm = TRUE)
bitsum = describeBy(mutresults$NeutralBits, mutresults$Experiment, mat=TRUE, na.rm = TRUE)
sizesum =  describeBy(mutresults$AvgGeneSize, mutresults$Experiment, mat=TRUE, na.rm = TRUE)
print(mutsum)
#dev.new()
#error.bars(stats=mutsum, labels=mutsum$group1, ylab="Neutral Mutations Rate" )

dev.new()
plot(mutsum$group1, mutsum$mean)
points(bitsum$group1, bitsum$mean)


print('###PROTEINS###')
sumpnum = describeBy(fresults$AvgNumProteins, fresults$Experiment, mat=TRUE)
sumfnum = describeBy(fresults$AvgNumFunctions, fresults$Experiment, mat=TRUE)
dev.new()
counts = data.frame(sumpnum$mean,sumfnum$mean)
print(counts)
barplot(t(counts), names.arg=sumpnum$group1,ylab="#Proteins/Functions", xlab='Experiment', main = "Number of Proteins and Functions in the Genomes", cex.names=CEXAXIS, col=c("darkblue","red"), legend = list('#Proteins','#Functions'), beside=TRUE, ylim=c(0,35))

dev.new()
plot(sumpnum$mean,sumresults$mean, ylab = "#Evaluations", xlab="#Proteins")

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
