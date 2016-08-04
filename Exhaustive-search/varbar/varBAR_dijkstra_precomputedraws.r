# varBAR_dijkstra_uniuni_precomputedraws.r

collect.draws<-function(pdflist) { #as function of indexer, aka pdflist, that which identifies the pdfs/distributions
    
    #set up matrix to store draws
    samples<-array(NA,dim=c(reps, length(pdflist)))
    colnames(samples)<-pdflist
    
    #set up progress checkpoints
    check01<-ceiling(length(pdflist)*0.01) #1% checkpoint
    check05<-ceiling(length(pdflist)*0.05) #5% checkpoint
    check10<-ceiling(length(pdflist)*0.1) #10% checkpoint
    checkcount<-0 #number of nodes sampled
    checkprint<-c(1,2,3,4,5,10,20,30,40,50,60,70,80,90,100)
    checkind<-1 #printing index
    
    for (i in 1:length(pdflist)) {
        print(paste('pre',i,sep=' '))
        #collect draws for the ith node
        drawcoords<-as.numeric(unlist(strsplit(pdflist[i],split='_')))
        drawpdf<-pdfmaker(drawcoords[1],drawcoords[2])
        #drawfunc<-function(x){drawpdf}
        drawfunc<-function(x){eval(parse(text=drawpdf[1]))}
        draws<-samplepdf(reps, drawfunc)
        samples[,i]<-draws
        print(paste('post',i,sep=' '))
        #write.table(samples, file=paste(filname,'-samplematrix.txt',sep=''), quote=F, sep='\t', row.names=F)
        
        #check progress
        checkcount<-checkcount+1
        if (checkcount<=check05 && checkcount%%check01==0) {
            print(date())
            print(paste('We are ', checkprint[checkind], '% done. States sampled:', sep=''))
            print(paste(i,'/',length(pdflist), sep=''))
            checkind<-checkind+1
        }
        if (checkcount%%check10==0) {
            print(date())
            print(paste('We are ', checkprint[checkind], '% done. States sampled:', sep=''))
            print(paste(i,'/',length(pdflist), sep=''))
            #write.table(samples, file=paste(filname,'-samplematrix.txt',sep=''), quote=F, sep='\t', row.names=F)
            #save(graph,file=paste(filname,'-map.RData',sep=''))
            checkind<-checkind+1
        }
    }
    
    # write.table(samples, file=paste(filname,'-samplematrix.txt',sep=''), quote=F, sep='\t', row.names=F)
    return(samples)
}
