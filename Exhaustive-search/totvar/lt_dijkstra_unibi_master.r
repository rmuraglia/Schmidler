# lt_dijkstra_unibi_master.r

# run dijkstra search for total variation distance in a grid indexed by lambda and temperature (here represented as sigma - sigma = sqrt(T))
# initial distn - normal
# target distn - bimodal

###########
# PART 1: set run parameters
###########

# store first time point for book keeping
time1<-proc.time()

# grab command line arguments if appropriate number are provided
args<-commandArgs(trailingOnly=TRUE)
if (length(args)==7) {
    filname<-as.character(args[1])
    move.jump<-as.numeric(args[2])
    reps<-as.numeric(args[3])
    sigma.min<-as.numeric(args[4])
    sigma.max<-as.numeric(args[5])
    sigma.numpoints<-as.numeric(args[6])
    lambda.numpoints<-as.numeric(args[7])
} else {
    filname<-'ltsearch-uni-bi-02'
    move.jump<-4
    reps<-1000
    sigma.min<-0.5
    sigma.max<-2.5
    sigma.numpoints<-9
    lambda.numpoints<-11
}


##########
# PART 2: print run info
##########

# create file sink for output. double output to console too
sink(file=paste(filname, '.sink', sep=''),split=T)

print(paste('Job started at ', date(), sep=''))
print(paste('Output files carry the ', filname, ' prefix.', sep=''))
print(paste('The move jump coefficient is ', move.jump, sep=''))
print(paste('Each distribution was sampled ', reps, ' times.', sep=''))

##########
# PART 3: run dijkstra's algorithm
##########

print('Importing auxillary scripts...')

source('lt_dijkstra_unibi_import.r')
source('lt_dijkstra_unibi_precomputedraws.r')

print(paste('Initial set up complete at ', date(), sep=''))
print(paste('number of nodes: ', nrow(point.grid), sep=''))
print(paste('number of potential moves per node: ', nrow(move.universe), sep=''))
print('initial state: ')
print(init)
print('target state: ')
print(target)
print('lambda values: ')
print(lambda.points)
print('sigma values: ')
print(sigma.points)

time2.1<-proc.time()

print('Performing sampling...')
samplematrix<-collect.draws(indexer)

time2.2<-proc.time()

print(paste('Sampling complete. Time for sampling was:'))
print(time2.2-time2.1)

print('Running Dijkstra\'s algorithm...')
ann.map<-do.dijkstra(map)

time2.3<-proc.time()
print(paste('Dijkstra\'s algorithm complete. Time for search was:'))
print(time2.3-time2.2)

#########
# PART 4: plot results
#########

print('Visualizing results')

# check solution path - it may need manual tuning before feeding to this final visualization
source('lt_dijkstra_unibi_viz.r')

time3<-proc.time()
print('Visualization complete. Time for viz was:')
print(time3-time2.3)

########
# PART 5: clean up and finish
########

time4<-proc.time()
print('R script complete')
print('Total time elapsed:')
print(time4-time1)

save.image(file=paste(filname, '.RData', sep=''))

print('Quitting R...')
sink()



