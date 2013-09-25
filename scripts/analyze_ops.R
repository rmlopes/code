library(psych)
library(xtable)
library(outliers)
library(tableplot)
library(ggplot2)
library(gridExtra)
options(width=192)
CEXAXIS <- 0.51

args<-commandArgs(TRUE)
separator = "\n\n"
#evoresults = read.table(args[1], header=TRUE)
evoresults = read.table("~/Documents/thesis/phdsupport/data/ops/analysis/grouped_evolution.txt", header=TRUE)
fresults = evoresults[evoresults$Problem != "harmonic" & (evoresults$Evaluations < 1000000 & evoresults$Problem != "<NA>" & !is.nan(evoresults$Best) & evoresults$Best < 0.001),]

#|| evoresults$Problem == 'harmonic',]
#fresults = fresults[fresults$Best < 0.01,]
#print(fresults[0:10])
fresults <- rbind(evoresults[evoresults$Problem == "harmonic",], fresults)
fresults$Problem <- factor(fresults$Problem)
fresults$Op <- factor(fresults$Op)
fresults$OpSize <- factor(fresults$OpSize)
fresults$Rates <- factor(fresults$Rates)
"
p1 <- ggplot(subset(fresults,Problem == 'harmonic'),
            aes(x = OpSize, y = Evaluations,fill=Rates)) + 
            facet_grid(Problem~Op) + 
            geom_boxplot()

p2 <- ggplot(subset(fresults,Problem == 'pendulum'),
            aes(x = OpSize, y = Evaluations,fill=Rates)) + 
            facet_grid(Problem~Op) + 
            geom_boxplot(outlier.shape=NA) +
            coord_cartesian(ylim=c(0, 50000))

p3 <- ggplot(subset(fresults,Problem == 'polinomial'),
            aes(x = OpSize, y = Evaluations,fill=Rates)) + 
            facet_grid(Problem~Op) + 
            geom_boxplot(outlier.shape=NA) +
            coord_cartesian(ylim=c(0, 300000))

p4 <- ggplot(subset(fresults,Problem == 'santafetrail'),
            aes(x=Problem,y = Evaluations,fill=Rates)) + 
            facet_grid(OpSize~Op) + 
            geom_boxplot(outlier.shape=NA) +
            coord_cartesian(ylim=c(0, 300000))
dev.new()
grid.arrange(p1,p2,p3,p4,ncol=1)
"

dev.new()
ggplot(data = fresults, aes(x = OpSize, y = Evaluations,fill=Rates)) + 
  geom_boxplot(outlier.shape=NA)+
  #coord_cartesian(ylim=c(0, 40000)) +
  facet_grid(Problem ~ Op ) +
  scale_y_log10() +
  theme(axis.text.x=element_text(angle=-90))


#splitresults = split(fresults,interaction(fresults$Problem,fresults$Op,fresults$Rates))
#print(splitresults)
#fresults$splitnames <- fresults$Problem.fresults$Op.fresults$Rates.fresults$OpSize
factors =  list(fresults$Problem,fresults$Op,fresults$Rates, fresults$OpSize)
sumresults = describeBy(fresults$Evaluations,interaction(fresults$Problem,fresults$Op,fresults$Rates, fresults$OpSize),mat=FALSE)
print(sumresults)
n <- factor(unique(interaction(fresults$Problem,fresults$Op,fresults$Rates, fresults$OpSize)))

#row.names(sumresults) <- levels(n)
#print(sumresults)
sumresultsb <- do.call("rbind",sumresults)
#row.names(sumresultsb) <- interaction(fresults$Problem,fresults$Op,fresults$Rates, fresults$OpSize)
sumresultsb$merged <- levels(n)
sumresultsb
#dev.new()
#b = boxplot(Evaluations ~interaction(Problem,Op), fresults,las=3, outline = FALSE, main="Number of Evaluations to succeed",col=(c("green","red",'blue','yellow')))
"ggplot(data = fresults, aes(x = factor(Rates), y = Evaluations, fill=Op)) + 
  geom_boxplot(outlier.shape=NA) +
  coord_cartesian(ylim=c(0, 500000)) +
  opts(axis.text.x=theme_text(angle=-90)) + facet_grid(Problem ~ OpSize) #+ 
  #opts(legend.position = 'none')
"
dev.new()
ggplot(data=fresults, aes(x=OpSize,fill=Rates)) +
  geom_bar(stat="bin",position="dodge") +
  #coord_flip() +
  facet_grid(Problem~Op) +
  theme(axis.text.x=element_text(angle=-90))

#factors = list(fresults$Rates,fresults$OpSize,fresults$Op, fresults$Problem)
#sumresults = describeBy(fresults$Evaluations, factors,na.rm = TRUE, mat=FALSE)
sumresults$Problem = strsplit(rowname)
dev.new()
error.bars(stats=sumresultsb, ylab="Evaluations" )

textable = xtable(sumresults[c(2,4:6,10)],
                           digits=c(0, 0, 0, 0,0,0),
				  caption="Summary of the results for the Symbolic Regression experiments.",
				  label="table:summary-symb")
print(textable, file = "summary-symb.tex")
print(sumresults)
cat(separator)
normcount = sumresults$n /100
dev.new()
#barplot(normcount, main='Successful Runs(%)', names.arg = sumresults$group1, cex.names=CEXAXIS)
ggplot(data = sumresults, aes(x = group1, y = n/100, fill=group2)) + 
  geom_bar()+ opts(axis.text.x=theme_text(angle=-90)) + facet_grid(group2 ~ .) + opts(legend.position = "none")

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
genresults = fresults[fresults$Validation > 0,]
genresults = genresults[genresults$Validation < 1000000,]
#genresults = fresults
sumgen =  describeBy(genresults$Validation, list(genresults$Experiment,genresults$Problem), mat=TRUE)
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
barplot(t(counts), names.arg=sumpnum$group1,ylab="#Proteins/Functions", xlab='Experiment', main = "Number of Proteins and Functions in the Genomes", las=3, col=c("darkblue","red"), legend = list('#Proteins','#Functions'), beside=TRUE, ylim=c(0,35))

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
