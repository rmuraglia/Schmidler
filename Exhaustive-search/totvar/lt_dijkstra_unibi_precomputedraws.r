# lt_dijkstra_unibi_precomputedraws.r

###########################
#PART 1: sampling engine
###########################

#functions for sampling taken from http://blog.quantitations.com/tutorial/2012/11/20/sampling-from-an-arbitrary-density/

endsign <- function(f, sign = 1) { #helper function to focus the search for uniroot
    b <- sign
    while (sign * f(b) < 0) { b <- b + sign*5 } #can tune how big to make search, depending on how spread the distribution is
    #while (sign * f(b) < 0) b <- 10 * b #original search expansion
    return(b)
}
#initially focus search in [-1,1], but if value isn't contained within, expand by adding +/- 5

samplepdf <- function(n, pdf, ..., spdf.lower = -Inf, spdf.upper = Inf) {
    vpdf <- function(v) sapply(v, pdf, ...)  # vectorize
    cdf <- function(z) integrate(vpdf, spdf.lower, z)$value
    invcdf <- function(u) {
        subcdf <- function(t) cdf(t) - u
        if (spdf.lower == -Inf) 
            spdf.lower <- endsign(subcdf, -1)
        if (spdf.upper == Inf) 
            spdf.upper <- endsign(subcdf)
        return(uniroot(subcdf, c(spdf.lower, spdf.upper))$root)
    }
    sapply(runif(n,min=0.01,max=0.99), invcdf)
}

#how does it work?
#1) samplepdf runs invcdf n times, with values randomly chosen from [0,1]
#2) in invcdf, endsign lines focus the search space
# invcdf takes in a random value from [0,1], and looks for what value of t inputted to the cdf is equal to that number.
# basically think of invcdf being given a value on the y axis. It projects it back onto the CDF curve and then projects down to find the value.
# the curvature/sharpness of the CDF determines the density along the x axis which is what we want.

#example call: 
#tim<-pdfmaker(0.5,1)
#samplepdf(1,function(x){eval(parse(text=tim[1]))})

##########################
#PART 2: distribution definitions
##########################

library(Ryacas)
Exp<-function(x) {exp(x)}
x<-Sym('x')

beta<-1
mu0<-0

U0<-function(k) { #define initial potential
    U <- (x-mu0)^2/(2*k^2)
    return(U)
}

U1<-function(k) { #define final potential
    U <- (x^4-81*x^2)/(150*k^2)
    return(U)
}

Umaker<-function(lambda,k) { #use lambda scaling to define intermediate potentials
    newU <- (1-lambda)*U0(k) + lambda*U1(k)
    return(newU)
}

Qget<-function(potential) { #calculate partition function. Usage: Qget(Umaker(0,1))
    integrand<-exp(-beta*potential)
    Q<-integrate(function(x){eval(parse(text=integrand[1]))},-Inf,Inf)
    return(Q)
}

pdfmaker<-function(lambda,k) {
    U<-Umaker(lambda,k)
    Q<-Qget(U)
    newpdf<-exp(-beta*U)/Q$value
    return(newpdf)
}

##########################
#PART 3: collect samples
##########################


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
    
    write.table(samples, file=paste(filname,'-samplematrix.txt',sep=''), quote=F, sep='\t', row.names=F)
    return(samples)
}