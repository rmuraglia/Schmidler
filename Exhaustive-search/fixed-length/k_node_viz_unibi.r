# k_node_viz_unibi.r

# import libraries
options(warn=-1)
suppressMessages(library(ggplot2))
suppressMessages(library(grid))
suppressMessages(library(gridExtra))
options(warn=0)

# set working directory
switch(Sys.info()[['sysname']],
    Windows = {setwd('D:/Dropbox/pathsearch-exhaustive/uni-bi-symmetric/')},
    Darwin = {setwd('/Users/rmuraglia/Dropbox/pathsearch-exhaustive/uni-bi-symmetric/')},
    Linux = {setwd('/dscrhome/rm204/pathsearch-exhaustive/uni-bi-symmetric/')}
)

if (!file.exists('imgout')) { dir.create('imgout') }

###########
# set up info panel
###########

emptydf<-data.frame(0,0)
colnames(emptydf)<-c('x', 'y')
labelplot<-ggplot(emptydf,aes(x=x,y=y)) + geom_blank() + theme_bw() + labs(x='',y='')
init_string<-paste('mu0=', mu0, ', sigma0=', sig0, sep='')
target_string<-paste('A=', a, ', B=', b, ', C=', c, sep='')
title_string<-paste('Best paths of fixed length \n', 'unimodal: ', init_string, '\nto bimodal: ', target_string, '\n distance metric: asympvarbar \nmove jump coeff: ', move_jump, '\nsamples per node: ', num_draws, sep='')

infopanel<-labelplot + geom_text(x=0, y=0, label=title_string, size=3) + theme(axis.ticks=element_blank(), axis.text.x=element_blank()) + theme(axis.ticks=element_blank(), axis.text.y=element_blank())

########
# individual path overviews
########

plotlist<-list(infopanel)
pcost_vec<-rep(NA, length.out=max_num_nodes)

for (i in 1:max_num_nodes) { # for each possible path length
    if (unscaled_cost[target_ind, i]!=Inf) { # if there exists a solution

        # get scaled path cost
        eff_N<-max_num_nodes*num_draws/i 
        path_cost<-unscaled_cost[target_ind, i]/(2*eff_N)
        print_cost<-format(signif(path_cost, digits=5), scientific=T)
        pcost_vec[i]<-path_cost

        # get path as data frame
        path_df<-path_soln(search_out[[3]], i)

        # make filename string
        i_print<-sprintf('%03g', i)
        imgname<-paste('imgout/', filID, '-', i_print, 'node.eps', sep='')

        # make plot
        baseA<-ggplot(data=point_grid, aes(x=lambda, y=temperature))
        layerA1<-geom_path(data=path_df, aes(x=lambda, y=temperature), arrow=arrow(), size=2, alpha=1)
        labelA1<-labs(title=paste(i, ' node path cost: ', print_cost, sep=''))
        plotA<-baseA + layerA1 + labelA1 + geom_point()
        plotlist[[i+1]]<-plotA

        # save to file
        ggsave(file=imgname, width=5, height=5, plot=plotA)
    }
}

##########
# summary figure for best paths
##########

top_ind<-order(pcost_vec, decreasing=F)[1:5]
summary_plots<-plotlist[c(1, top_ind+1)]
plotB<-arrangeGrob(summary_plots[[1]], summary_plots[[2]], summary_plots[[3]], summary_plots[[4]], summary_plots[[5]], summary_plots[[6]], nrow=2)
summaryname<-paste('imgout/', filID, '-summary.eps', sep='')
ggsave(file=summaryname, width=9, height=6, plot=plotB)