# varBAR_dijkstra_uniuni_import.r

# import dijkstra methods

##########################
#PART 1: Grid setup
##########################

#specify ranges and values for temperature and lambda
if (!exists('temp.min')) {temp.min<-0.5} #0.5, 1
if (!exists('temp.max')) {temp.max<-2.5} #2.5, 2
if (!exists('temp.numpoints')) {temp.numpoints<-9} #9, 3
temp.points<-seq(from=temp.min, to=temp.max, length.out=temp.numpoints)

lambda.min<-0
lambda.max<-1
if (!exists('lambda.numpoints')) {lambda.numpoints<-11} #11, 6
lambda.points<-seq(from=lambda.min, to=lambda.max, length.out=lambda.numpoints)

#select inital and target coordinates
init<-c(lambda.points[1],temp.points[3])
target<-c(lambda.points[length(lambda.points)],temp.points[3])

#create objects holding point set information
point.lists<-list(lambda.points,temp.points)
names(point.lists)<-c('lambda','temp')
point.grid<-expand.grid(point.lists)
point.list<-as.list(as.data.frame(t(point.grid)))

########################
#PART 2: Dijkstra map setup
########################

#create map to be annotated during Dijkstra search
map<-list()
for (i in 1:length(point.list)) {
    map[[length(map)+1]]<-list(coords=point.list[[i]], totdist=1E10, visit=FALSE, prev=NA)#, samples=NA)
}

#create indexing tool for map navigation
indexer<-apply(point.grid,1,paste,collapse='_')

init.key<-paste(init,collapse='_')
init.ind<-which(indexer==init.key)

#initialize initial state cumulative distance to zero
map[[init.ind]]$totdist<-0

########################
#PART 3: Movement definitions
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

#movement primitives
move.result<-function(move,state) {
    #compute new state based on original state and set of changes from the move
    new.pos<-rep(NA,length.out=length(state))
    
    for (i in 1:length(move)) { #for each dimension
        if (!is.na(move[i])) { #if free to move
        #get potential values
            dimension<-point.lists[[i]]
            curr.ind<-which(dimension==state[i]) #find where it is
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

########################
#PART 4: Dijkstra search routines
########################

do.dijkstra<-function(graph) {
    
    #set up checkpoints for run-time diagnostics
    check01<-ceiling(length(graph)*0.01) #1% checkpoint
    check05<-ceiling(length(graph)*0.05) #5% checkpoint
    check10<-ceiling(length(graph)*0.1) #10% checkpoint
    checkcount<-0 #count number of popped/visited nodes as measure of progress
    checkprint<-c(1,2,3,4,5,10,20,30,40,50,60,70,80,90,100)
    checkind<-1
    
    #while there is at least one node that has not been visited, continue searching
    while(any(!sapply(graph,'[[',3))) {
        
        ###############
        #We want to pop off the node that has the lowest best distance that has not yet been visited, and then expand over their neighbors
        ###############
        
        #get list of nodes that have already been visited
        visited<-which(sapply(graph,'[[',3))
        
        #get current best distances
        qdist<-sapply(graph,'[[',2)
        
        #if already visited, we do not want to visit it again: overwrite current best distance with something nonsensical
        qdist[visited]<-NA
        
        #get index of the lowest distance in the queue
        mindist<-which.min(qdist)
        
        #get node information for that min distance node
        u<-graph[[mindist]]
        
        #marking it as visited it equivalent to popping it off the queue: it will not be visited again
        graph[[mindist]]$visit<-T
        
        #check our progress
        checkcount<-checkcount+1 #update counter
        #print(checkcount)
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
            #save(graph,file=paste(filname,'-map.RData',sep=''))
            checkind<-checkind+1
        }
        
        ##################
        #now that we have a node to inspect, expand out over its neighbors, and see if we can improve paths to them
        ##################
        
        #get neighbors of u
        move.bool<-apply(move.universe, 1, is.valid.move, u$coords)
        valid.moves<-move.universe[move.bool,]
        neighbors<-t(apply(valid.moves, 1, move.result, u$coords))
        
        #for each neighbor, calculate total distance to that neighbor, going through u
        for (i in 1:nrow(neighbors)) {
            
            neighbor<-neighbors[i,] #get neighbor
            neighdist<-dist.move(u$coords, neighbor) #get dist from u to neighbor
            alt<-u$totdist + neighdist #total distance to neighbor through u
            
            #compare alternate distance to current best for the neighbor
            v.ind<-which(indexer==paste(neighbor, collapse='_'))
            v.dist<-graph[[v.ind]]$totdist
            
            #if alt distance is better, then update this node's info
            if (alt < v.dist) {
                graph[[v.ind]]$totdist<-alt
                graph[[v.ind]]$prev<-u$coords
            }
        } #close loop over neighbors of u
    } #close loop over queue
    
    return(graph)
}

dijkstra.soln<-function(graph) {
    
    #initialize the output sequence as the target coordinates
    seq.out<-list(target)
    
    #get the index of that node in the graph list
    trace.index<-which(indexer==paste(target,collapse='_'))
    
    #keep tracking back until the listed previous node is an NA (meaning it has no predecessor)
    while (!any(is.na(graph[[trace.index]]$prev))) {
        
        #add the coordinates of the previous node to the beginning of the output sequence
        seq.out<-append(seq.out,list(graph[[trace.index]]$prev),after=0)
        
        #get the next previous node
        trace.index<-which(indexer==paste(graph[[trace.index]]$prev,collapse='_'))
    }
    
    #return a list of coordinates, in sequence from init to target
    return(seq.out)
}


path.cost<-function(path) { #get total cost of a path
    cost<-0
    for (i in 2:length(path)) {
        A<-path[[i-1]]
        B<-path[[i]]
        cost<-cost+dist.move(A,B)
    }
    return(cost)
}

totalratio<-function(path) { #get product of ratios along a path
    ratios<-rep(NA,length.out=(length(path)-1))
    for (i in 2:length(path)) {
        A<-path[[i-1]]
        B<-path[[i]]
        ratio<-ratio.estim(A,B,samplematrix)
        ratios[i-1]<-ratio
    }
    return(ratios)
}