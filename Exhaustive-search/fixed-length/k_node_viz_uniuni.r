# k_node_viz_uniuni.r

###########
# info panel
###########

emptydf<-data.frame(0,0)
colnames(emptydf)<-c('x','y')
labelplot<-ggplot(emptydf,aes(x=x,y=y)) + geom_blank() + theme_bw() + labs(x='',y='')
init_string<-paste('mu0=', mu0, ', sigma0=', sig0, sep='')
target_string<-paste('mu1=', mu1, ', sigma1=', sig1, sep='')
title_string<-paste('Best paths of fixed length \n', '1D normal \n', init_string, '\n', target_string, '\n Distance metric: asympvarbar \n Move jump coeff: ', move.jump, '\n Samples per node: ', reps, sep='')

infopanel<-labelplot + geom_text(x=0 ,y=0, label=title_string, size=3) + theme(axis.ticks = element_blank(), axis.text.x = element_blank()) +theme(axis.ticks = element_blank(), axis.text.y = element_blank())

#########
# individual path overviews
#########
plotlist<-list(infopanel)
pcost_vec<-rep(NA, length.out=max_num_nodes)

for (i in max_num_nodes:1) {
    if (unscaled_cost[target.ind, i]!=Inf) {

        # get scaled path cost
        eff_N<-nmax/i
        path_cost<-unscaled_cost[target.ind, i]/(2*eff_N)
        print_cost<-format(signif(path_cost, digits=5), scientific=T)
        pcost_vec[i]<-path_cost

        # get path as data frame
        path_df<-path_soln(search_out[[3]], i)

        # make filename string
        imgname<-paste(filname, '-', i, 'node.png', sep='')

        # make plot
        p1<-ggplot(data=point.grid, aes(x=lambda, y=temp))
        a1<-geom_path(data=path_df, aes(x=lambda, y=temp), arrow=arrow(), size=2, alpha=1)
        a2<-labs(title=paste(i, ' node path cost: ', print_cost, sep=''))
        full_plot<-p1+a1+a2+geom_point()
        plotlist[[i+1]]<-full_plot

        # save plot
        png(imgname, width=3, height=3, units='in', res=200)
        print(full_plot)
        dev.off()
    }
}

#####
# collective summary of top 5 
#####

top_ind<-order(pcost_vec, decreasing=F)[1:5]
summary_plots<-plotlist[c(1, top_ind+1)]
png(paste(filname, '-summary.png', sep=''), width=10, height=6, units='in', res=200)
do.call(grid.arrange, c(summary_plots, list(nrow=2)))
dev.off()