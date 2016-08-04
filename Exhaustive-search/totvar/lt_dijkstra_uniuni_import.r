# lt_dijkstra_uniuni_import.r

# import functions and build objects for dijkstra search
# total variation as distance metric, unimodal end distributions.

##########################
#PART 1: Grid setup
##########################

# specify ranges and values for sigma and lambda
if (!exists('sigma.min')) { sigma.min<-0.5 } # set default if master fails to define
if (!exists('sigma.max')) { sigma.max<-2.5 } # set default if master fails to define
if (!exists('sigma.numpoints')) { sigma.numpoints<-9 } # set default if master fails to define
sigma.points<-seq(from=sigma.min, to=sigma.max, length.out=sigma.numpoints)

lambda.min<-0
lambda.max<-1
if (!exists('lambda.numpoints')) { lambda.numpoints<-11 } # set default if master fails to define
lambda.points<-seq(from=lambda.min, to=lambda.max, length.out=lambda.numpoints)

# select inital and target coordinates
init<-c(lambda.points[1], sigma.points[3])
target<-c(lambda.points[length(lambda.points)], sigma.points[3])

# create objects holding point set information
point.lists<-list(lambda.points, sigma.points)
names(point.lists)<-c('lambda', 'sigma')
point.grid<-expand.grid(point.lists)
point.list<-as.list(as.data.frame(t(point.grid)))

########################
#PART 2: Dijkstra map setup
########################

# create map to be annotated during Dijkstra search
map<-list()
for (i in 1:length(point.list)) {
    map[[length(map)+1]]<-list(coords=point.list[[i]], totdist=1E10, visit=FALSE, prev=NA, samples=NA)
}

# create indexing tool for map navigation
indexer<-apply(point.grid, 1, paste, collapse='_')

init.key<-paste(init, collapse='_')
init.ind<-which(indexer==init.key)

# initialize initial state cumulative distance to zero
map[[init.ind]]$totdist<-0

########################
#PART 3: Movement definitions
########################

# movement options
dof<-c(TRUE,TRUE) # allow both lambda and sigma to vary - if want restricted, fixed sigma path, set dof<-c(TRUE, FALSE)
names(dof)<-c('lambda','sigma')

dof.ranges<-list()
for (i in dof) {
    if (i) { dof.ranges<-append(dof.ranges,list(-move.jump:move.jump))
    } else { dof.ranges<-append(dof.ranges,list(NA)) }
}
names(dof.ranges)<-names(dof)

nacheck<-function(x) { return(any(is.na(x))) }

get.all.moves<-function() {
    grid.full<-expand.grid(dof.ranges)
    grid.strip<-grid.full[,!apply(grid.full,2,nacheck)] # remove NA columns that mess with absolute value
    grid.trim1<-grid.strip[rowSums(abs(grid.strip))<=move.jump,] # filter out actions with too many moves
    grid.trim2<-grid.trim1[rowSums(abs(grid.trim1))>0,] # filter out null moves
    ind<-as.numeric(rownames(grid.trim2))
    grid.out<-grid.full[ind,]
    return(grid.out)
}

move.universe<-get.all.moves()

# movement primitives
move.result<-function(move,state) {
    # compute new state based on original state and set of changes from the move
    new.pos<-rep(NA,length.out=length(state))
    
    for (i in 1:length(move)) { # for each dimension
        # get potential values
        dimension<-point.lists[[i]]
        curr.ind<-which(dimension==state[i]) # find where it is
        new.ind<-curr.ind + as.numeric(move[i]) # get new index
        if (new.ind<=0) { new.pos[i]<-NA # if out of bounds, call it NA
        } else if (new.ind>length(dimension)) { new.pos[i]<-NA
        } else {new.pos[i]<-dimension[new.ind] } # otherwise get new position
    } # note: output may have NAs in position: this is fine, since it will get pruned by the valid move checker
    return(new.pos)
}

is.valid.move<-function(move,state) {
    pos<-paste(move.result(move,state),collapse='_')
    if (pos %in% indexer) { return(TRUE)
    } else { return(FALSE) }
}

########################
#PART 4: Dijkstra search routines
########################

do.dijkstra<-function(graph) {
    
    # set up checkpoints for run-time diagnostics
    check01<-ceiling(length(graph)*0.01) # 1% checkpoint
    check05<-ceiling(length(graph)*0.05) # 5% checkpoint
    check10<-ceiling(length(graph)*0.10) # 10% checkpoint
    checkcount<-0 # count number of popped/visited nodes as measure of progress
    checkprint<-c(1,2,3,4,5,10,20,30,40,50,60,70,80,90,100)
    checkind<-1
    
    # while there is at least one node that has not been visited, continue searching
    while(any(!sapply(graph,'[[',3))) {
        
        ###############
        # We want to pop off the node that has the lowest best distance that has not yet been visited, and then expand over their neighbors
        ###############
        
        # get list of nodes that have already been visited
        visited<-which(sapply(graph,'[[',3))
        
        # get current best distances
        qdist<-sapply(graph,'[[',2)
        
        # if already visited, we do not want to visit it again: overwrite current best distance with something nonsensical
        qdist[visited]<-NA
        
        # get index of the lowest distance in the queue
        mindist<-which.min(qdist)
        
        # get node information for that min distance node
        u<-graph[[mindist]]
        
        # marking it as visited it equivalent to popping it off the queue: it will not be visited again
        graph[[mindist]]$visit<-T
        
        ######
        # draw samples for this node: two options
        ######

        # option 1: general purpose for all pdfs
        updf<-pdfmaker(u$coords[1],u$coords[2])
        updffunc<-function(x){eval(parse(text=updf[1]))}
        udraws<-samplepdf(reps,updffunc)
        u$samples<-udraws
        graph[[mindist]]$samples<-udraws

        # option 2: specific to normal-normal example
        # all distributions in grid will be normal, so derive updated params and use normal sampler
        # note: this version may be slightly off due to how distributions are defined in this old scheme
        # note: likely requires an alternate totvar dist calculator which uses pnorm with the param updates
        # mu_new<-mu_update(u$coords[1], u$coords[2]^2)
        # sig_new<-sig_update(u$coords[1], u$coords[2]^2)
        # udraws<-norm_sampler(mu_new, sig_new, reps)
        # u$samples<-udraws
        # graph[[mindist]]$samples<-udraws

        #####
        # end draw samples
        #####
        
        # check our progress
        checkcount<-checkcount+1 # update counter
        print(checkcount)
        if (checkcount<=check05 && checkcount%%check01==0) {
            print(date())
            print(paste('We are ', checkprint[checkind], '% done. Popped node progress:', sep=''))
            print(table(sapply(graph,'[[',3)))
            checkind<-checkind+1
        }
        if (checkcount%%check10==0) {
            print(date())
            print(paste('We are ', checkprint[checkind], '% done. Popped node progress:', sep=''))
            print(table(sapply(graph,'[[',3)))
            save(graph,file=paste(filname,'-map.RData',sep=''))
            checkind<-checkind+1
        }
        
        ##################
        # now that we have a node to inspect, expand out over its neighbors, and see if we can improve paths to them
        ##################
        
        # get neighbors of u
        move.bool<-apply(move.universe, 1, is.valid.move, u$coords)
        valid.moves<-move.universe[move.bool,]
        neighbors<-t(apply(valid.moves, 1, move.result, u$coords))
        
        # for each neighbor, calculate total distance to that neighbor, going through u
        for (i in 1:nrow(neighbors)) {
            
            neighbor<-neighbors[i,] # get neighbor
            neighdist<-totvar.dist(u, neighbor) # get dist from u to neighbor
            alt<-u$totdist + neighdist # total distance to neighbor through u
            
            # compare alternate distance to current best for the neighbor
            v.ind<-which(indexer==paste(neighbor, collapse='_'))
            v.dist<-graph[[v.ind]]$totdist
            
            # if alt distance is better, then update this node's info
            if (alt < v.dist) {
                graph[[v.ind]]$totdist<-alt
                graph[[v.ind]]$prev<-u$coords
            }
        } # close loop over neighbors of u
    } # close loop over queue
    
    return(graph)
}

# take annotated map and get solution path

dijkstra.soln<-function(graph) {
    
    # initialize the output sequence as the target coordinates
    seq.out<-list(target)
    
    # get the index of that node in the graph list
    trace.index<-which(indexer==paste(target,collapse='_'))
    
    # keep tracking back until the listed previous node is an NA (meaning it has no predecessor)
    while (!any(is.na(graph[[trace.index]]$prev))) {
        
        # add the coordinates of the previous node to the beginning of the output sequence
        seq.out<-append(seq.out,list(graph[[trace.index]]$prev),after=0)
        
        # get the next previous node
        trace.index<-which(indexer==paste(graph[[trace.index]]$prev,collapse='_'))
    }
    
    # return a list of coordinates, in sequence from init to target
    return(seq.out)
}

########################
# PART 5: Endpoint distribution definitions and helper functions
########################

library(Ryacas)
Exp<-function(x) {exp(x)}
x<-Sym('x')

beta<-1
mu0<-0
mu1<-5

U0<-function(k) { # define initial potential
    U <- (x-mu0)^2/(2*k^2)
    return(U)
}

U1<-function(k) { # define final potential
    U <- (x-mu1)^2/(2*k^2)
    return(U)
}

Umaker<-function(lambda,k) { # use lambda scaling to define intermediate potentials
    newU <- (1-lambda)*U0(k) + lambda*U1(k)
    return(newU)
}

Qget<-function(potential) { # calculate partition function. Usage: Qget(Umaker(0,1))
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

###########################
# PART 6: Sampling method
###########################

# option 1: general for any pdf

# functions for sampling taken from http://blog.quantitations.com/tutorial/2012/11/20/sampling-from-an-arbitrary-density/

endsign <- function(f, sign = 1) { # helper function to focus the search for uniroot
    b <- sign
    while (sign * f(b) < 0) { b <- b + sign*5 } # can tune how big to make search, depending on how spread the distribution is
    #while (sign * f(b) < 0) b <- 10 * b #original search expansion
    return(b)
}
# initially focus search in [-1,1], but if value isn't contained within, expand by adding +/- 5

samplepdf <- function(n, pdf, ..., spdf.lower = -Inf, spdf.upper = Inf) {
    vpdf <- function(v) sapply(v, pdf, ...)  # vectorize
    cdf <- function(x) integrate(vpdf, spdf.lower, x)$value
    invcdf <- function(u) {
        subcdf <- function(t) cdf(t) - u
        if (spdf.lower == -Inf) 
            spdf.lower <- endsign(subcdf, -1)
        if (spdf.upper == Inf) 
            spdf.upper <- endsign(subcdf)
        return(uniroot(subcdf, c(spdf.lower, spdf.upper))$root)
    }
    sapply(runif(n), invcdf)
}

# how does it work?
# 1) samplepdf runs invcdf n times, with values randomly chosen from [0,1]
# 2) in invcdf, endsign lines focus the search space
# invcdf takes in a random value from [0,1], and looks for what value of t inputted to the cdf is equal to that number.
# basically think of invcdf being given a value on the y axis. It projects it back onto the CDF curve and then projects down to find the value.
# the curvature/sharpness of the CDF determines the density along the x axis which is what we want.

# example call: 
# tim<-pdfmaker(0.5,1)
# samplepdf(1,function(x){eval(parse(text=tim[1]))})

# option 2: specific to normals example
# sig1<-target[2]
# sig0<-init[2]

# mu_update<-function(lambda, T) {
#     mu_new<-(sig1^2*mu0*(1-lambda) + sig0^2*mu1*lambda) / (sig1^2*(1-lambda) + sig0^2*lambda)
#     return(mu_new)
# }

# # sigma update rule for normals:
# sig_update<-function(lambda, T) {
#     sig_new<-sqrt( T*(sig0^2*sig1^2) / (sig1^2*(1-lambda) + sig0^2*lambda) )
#     return(sig_new)
# }

# norm_sampler<-function(mu, sigma, ndraws) { rnorm(ndraws, mu, sigma) }

################################
#PART 7: distance evaluation
################################

# source for this expression of total variation distance:
# Gareth O Roberts and Jeffrey S Rosenthal. General state space Markov chains and MCMC algorithms. Probability Surveys, 1:20–71, 2004.
# d_{TV}(p_1, p_2) = E_1[min(1, p_2(x)/p_1(x))]

totvar.dist<-function(node,state2) {
    state1<-node$coords
    
    pdf1<-pdfmaker(state1[1],state1[2])
    pdf2<-pdfmaker(state2[1],state2[2])
    pdf1func<-function(x){eval(parse(text=pdf1[1]))}
    pdf2func<-function(x){eval(parse(text=pdf2[1]))}
    
    draws<-node$samples
    dens1<-pdf1func(draws)
    dens2<-pdf2func(draws)
    
    a<-dens2/dens1
    b<-which(a>1)
    a[b]<-1
    repout<-1-a
    return(mean(repout))
}