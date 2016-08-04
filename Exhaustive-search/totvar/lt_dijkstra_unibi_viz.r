# lt_dijkstra_unibi_viz.r

library(ggplot2)
library(grid)
library(gridExtra)
library(gplots)

######################
#PART 1: Calculate path info
######################

soln.path<-dijkstra.soln(ann.map)

path.cost<-function(path) { #get total cost of a path
    cost<-0
    for (i in 2:length(path)) {
        A<-path[[i-1]]
        B<-path[[i]]
        cost<-cost+totvar.dist(A,B)
    }
    return(cost)
}

dcost<-path.cost(soln.path)
dpath<-as.data.frame(do.call(rbind,soln.path))
colnames(dpath)<-c('lambda','sigma')
soln.indexer<-apply(dpath,1,paste,collapse='_')

########
# PART 2: set up optimal path density plots
########

den.lim<-0.01 #density cutoff for evaluation
lim.l<--Inf
lim.u<-Inf

cdf2<-function(x, pdf) {integrate(pdf,lim.l,x)$value }
invcdf2<-function(u, pdf) {
    subcdf<-function(t) {cdf2(t, pdf)-u}
    if (lim.l == -Inf) { lim.l<-endsign(subcdf,-1) }
    if (lim.u == Inf) { lim.u<-endsign(subcdf) }
    return(uniroot(subcdf, c(lim.l,lim.u))$root)
}

lowers<-rep(NA,length.out=length(soln.indexer))
uppers<-lowers

for (i in 1:length(soln.indexer)) {
    A<-as.numeric(unlist(strsplit(soln.indexer[i],split='_')))
    pdfA<-pdfmaker(A[1],A[2])
    funcA<-function(x) {eval(parse(text=pdfA[1]))}
    lowers[i]<-invcdf2(den.lim,funcA)
    uppers[i]<-invcdf2(1-den.lim, funcA)
}

xvals<-seq(from=min(lowers), to=max(uppers), length.out=500)

plotdat<-array(NA, dim=c((length(xvals)*length(soln.indexer)),3))
colnames(plotdat)<-c('values','density','ind')

for (i in 1:length(soln.indexer)) {
    A<-as.numeric(unlist(strsplit(soln.indexer[i],split='_')))
    pdfA<-pdfmaker(A[1],A[2])
    funcA<-function(x) {eval(parse(text=pdfA[1]))}
    Avals<-funcA(xvals)
    plotdat[(1+(i-1)*length(xvals)):(i*length(xvals)),1]<-xvals
    plotdat[(1+(i-1)*length(xvals)):(i*length(xvals)),2]<-Avals
    plotdat[(1+(i-1)*length(xvals)):(i*length(xvals)),3]<-soln.indexer[i]
}

plot.df<-as.data.frame(plotdat)
plot.df[,1]<-as.numeric(as.character(plot.df[,1]))
plot.df[,2]<-as.numeric(as.character(plot.df[,2]))
plot.df[,3]<-factor(plot.df[,3],levels=soln.indexer) #put all factors in order of path

pcolors<-rich.colors(length(soln.indexer))
p1<-ggplot(data=plot.df, aes(x=values, y=density))
a1<-geom_area(aes(colour=ind, fill=ind), alpha=0.3, position='identity',size=1)
a2<-scale_fill_manual(values=pcolors, name='Optimal\ndistribution\nsequence')
a3<-scale_colour_manual(values=pcolors, name='Optimal\ndistribution\nsequence')
a4<-labs(title=paste('Augmented lambda-scaling path cost: ', round(path.cost(soln.path), 3)))

####### 
# PART 3: set up standard (bad) path viz
#######

bad.indexer<-c('0_1','0.1_1','0.2_1','0.4_1','0.5_1','0.6_1','0.8_1','0.9_1','1_1')

baddat<-array(NA, dim=c((length(xvals)*length(bad.indexer)),3))
colnames(baddat)<-c('values','density','ind')

for (i in 1:length(bad.indexer)) {
    A<-as.numeric(unlist(strsplit(bad.indexer[i],split='_')))
    pdfA<-pdfmaker(A[1],A[2])
    funcA<-function(x) {eval(parse(text=pdfA[1]))}
    Avals<-funcA(xvals)
    baddat[(1+(i-1)*length(xvals)):(i*length(xvals)),1]<-xvals
    baddat[(1+(i-1)*length(xvals)):(i*length(xvals)),2]<-Avals
    baddat[(1+(i-1)*length(xvals)):(i*length(xvals)),3]<-bad.indexer[i]
}

bad.df<-as.data.frame(baddat)
bad.df[,1]<-as.numeric(as.character(bad.df[,1]))
bad.df[,2]<-as.numeric(as.character(bad.df[,2]))
bad.df[,3]<-factor(bad.df[,3],levels=bad.indexer) #put all factors in order of path

bad.path<-list(c(0,1)); bad.path[[2]]<-c(0.1,1); bad.path[[3]]<-c(0.2,1)
bad.path[[4]]<-c(0.4,1); bad.path[[5]]<-c(0.5,1); bad.path[[6]]<-c(0.6,1)
bad.path[[7]]<-c(0.8,1); bad.path[[8]]<-c(0.9,1); bad.path[[9]]<-c(1,1)

b1<-ggplot(data=bad.df, aes(x=values, y=density))
d1<-geom_area(aes(colour=ind, fill=ind), alpha=0.3, position='identity',size=1)
d2<-scale_fill_manual(values=pcolors, name='Suboptimal\ndistribution\nsequence')
d3<-scale_colour_manual(values=pcolors, name='Suboptimal\ndistribution\nsequence')
d4<-labs(title=paste('Standard lambda-scaling path cost: ', round(path.cost(bad.path), 3)))

##########
# PART 4: parameter overview plot
##########

dpath<-cbind(dpath,'Augmented')
colnames(dpath)[3]<-'ind'
bpath<-dpath
bpath[,1]<-c(0,0.1,0.2,0.4,0.5,0.6,0.8,0.9,1)
bpath[,2]<-1
bpath[,3]<-'Standard'
path.df<-rbind(dpath,bpath)
pp1<-geom_path(data=path.df, aes(x=lambda, y=sigma, colour=ind), arrow=arrow(), size=2)
pp2<-scale_colour_discrete(name=' ')
pp3<-theme(legend.position='bottom')
paramplot2<-ggplot(data=point.grid, aes(x=lambda, y=sigma)) + pp1 + pp2 + pp3 + geom_point()

#########
# PART 5: put it all together and output
#########

png('unib1-3pan-searchv1.png', width=15, height=6, units='in',res=600)
grid.arrange(paramplot2+coord_fixed(ratio=0.5), arrangeGrob(p1+a1+a2+a3+a4, b1+d1+d2+d3+d4, nrow=2), ncol=2, widths=c(2,5))
dev.off()
