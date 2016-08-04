# grid_import_uniuni.r

##########################
#PART 1: Grid setup
##########################

# set temperature points
if (!exists('temp.min')) {temp.min<-0.5} 
if (!exists('temp.max')) {temp.max<-2.5} 
if (!exists('temp.numpoints')) {temp.numpoints<-9} 
temp.points<-seq(from=temp.min, to=temp.max, length.out=temp.numpoints)
temp.ind<-which.min(abs(temp.points-1.5))


# set lambda point values
lambda.min<-0
lambda.max<-1
if (!exists('lambda.numpoints')) { lambda.numpoints<-11 }
lambda.points<-seq(from=lambda.min, to=lambda.max, length.out=lambda.numpoints)

# set initial and target states
init<-c(lambda.points[1],temp.points[temp.ind])
target<-c(lambda.points[lambda.numpoints],temp.points[temp.ind])

# create objectes holding point set information
point.lists<-list(lambda.points, temp.points)
names(point.lists)<-c('lambda', 'temp')
point.grid<-expand.grid(point.lists)
point.list<-as.list(as.data.frame(t(point.grid)))
indexer<-apply(point.grid,1,paste,collapse='_')
init.key<-paste(init,collapse='_')
init.ind<-which(indexer==init.key)
target.key<-paste(target, collapse='_')
target.ind<-which(indexer==target.key)

########################
#PART 2: Movement definitions
########################

#movement options
# dof<-c(TRUE,TRUE)
# dof<-c(TRUE,FALSE)
names(dof)<-c('lambda','temp')

dof.ranges<-list()
for (i in dof) {
    if (i) { dof.ranges<-append(dof.ranges,list(-move.jump:move.jump))
    } else { dof.ranges<-append(dof.ranges,list(NA)) }
}
names(dof.ranges)<-names(dof)

nacheck<-function(x) { return(any(is.na(x))) }

get.all.moves<-function() {
    grid.full<-expand.grid(dof.ranges)
    grid.strip<-grid.full[,!apply(grid.full,2,nacheck),drop=F] #remove NA columns that mess with absolute value
    grid.trim1<-grid.strip[rowSums(abs(grid.strip))<=move.jump,,drop=F] #filter out actions with too many moves
    grid.trim2<-grid.trim1[rowSums(abs(grid.trim1))>0,,drop=F] #filter out null moves
    ind<-as.numeric(rownames(grid.trim2))
    grid.out<-grid.full[ind,]
    return(grid.out)
}

move.universe<-get.all.moves()

# helper function to replace numeric comparison with ==
equivcheck<-function(list, point) { isTRUE(all.equal(list, point)) }

#movement primitives
move.result<-function(move,state) {
    #compute new state based on original state and set of changes from the move
    new.pos<-rep(NA,length.out=length(state))
    
    for (i in 1:length(move)) { #for each dimension
        if (!is.na(move[i])) { #if free to move
        #get potential values
            dimension<-point.lists[[i]]
            curr.ind<-which(sapply(dimension, equivcheck, state[i]))
            # curr.ind<-which(dimension==state[i]) #find where it is #legacy shitty version that used == to compare numerics
            new.ind<-curr.ind + as.numeric(move[i]) #get new index
            if (new.ind<=0) { new.pos[i]<-NA #if out of bounds, call it NA
            } else if (new.ind>length(dimension)) { new.pos[i]<-NA
            } else {new.pos[i]<-dimension[new.ind] } #otherwise get new position
        } else { new.pos[i]<-state[i] } #if move is an NA then just stay where you are
    } #note: output may have NAs in position: this is fine, since it will get pruned by the valid move checker
    return(new.pos)
}

is.valid.move<-function(move,state) {
    pos<-paste(move.result(move,state),collapse='_')
    if (pos %in% indexer) { return(TRUE)
    } else { return(FALSE) }
}

####################
# PART 3: Sample collection routine
####################

collect_draws<-function(state, numsam) {
    lambda<-state[1]
    temp<-state[2]

    #get parameter updates
    mu_new<-mu_update(lambda, temp)
    sig_new<-sig_update(lambda, temp)

    #draw samples
    samps<-rnorm(n=numsam, mean=mu_new, sd=sig_new)
    return(samps)
}