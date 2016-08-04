# k_node_unibi_master.r

# run path search on varBAR distance in a grid indexed by lambda and temperature
# initial distribution is normal
# target distribution is bimodal
# here we find paths of fixed length, k
# this allows us to adjust the path cost to account for the total number of edges in the path, which is a proxy of computational cost of traversing that path.

# solution based on http://stackoverflow.com/questions/1690347/shortest-path-with-a-fixed-number-of-edges

###############
# PART 1: Set up and initialization
###############

print('begin R script')
print(date())

print('reading parameters from command line')

# set working directory
switch(Sys.info()[['sysname']],
    Windows = {setwd('D:/Dropbox/pathsearch-exhaustive/uni-bi-symmetric/')},
    Darwin = {setwd('/Users/rmuraglia/Dropbox/pathsearch-exhaustive/uni-bi-symmetric/')},
    Linux = {setwd('/dscrhome/rm204/pathsearch-exhaustive/uni-bi-symmetric/')}
)

# grab command line args if there are any
args<-commandArgs(trailingOnly=TRUE)

if (length(args)==7) { # if appropriate number, use
    move_jump<-as.numeric(args[1])
    temp_numpoints<-as.numeric(args[2])
    lambda_numpoints<-as.numeric(args[3])
    temp_min<-as.numeric(args[4])
    temp_max<-as.numeric(args[5])
    temp_dof<-as.logical(args[6])
    num_draws<-as.numeric(args[7])
} else { # else use defaults
    move_jump<-10
    temp_numpoints<-10
    lambda_numpoints<-11
    temp_min<-0.5
    temp_max<-5
    temp_dof<-TRUE
    num_draws<-1000
}

# create identifiers
filID<-paste('pathsearch-', temp_numpoints, 'tpts-', lambda_numpoints, 'lpts-', move_jump, 'jump-', temp_dof, 'tdof-', num_draws, '-draws', sep='')
gridID<-paste(temp_numpoints,'tpts-', lambda_numpoints, 'lpts-', temp_min, 'tmin-', temp_max, 'tmax', sep='')

print(filID)

# if exists, load precomputed distribution objects from ryacas
if (file.exists(paste('yacas-import-', gridID, '.RData', sep=''))) {
    load(paste('yacas-import-', gridID, '.RData', sep=''))
} else { 
    print('ERROR: no precomputed yacas object found') 
    print('please run yacas_precompute.r to create a yacas object for your grid')
    quit('no')
}

# if samplematrix exists, load it up, else sample it
if (file.exists(paste('samplematrix-', gridID, '-', num_draws, 'numdraws', '.txt', sep=''))) {
    print('loading precomputed sample matrix:')
    print(paste('samplematrix-', gridID, '-', num_draws, 'numdraws', '.txt', sep=''))
    sample_matrix<-read.table(paste('samplematrix-', gridID, '-', num_draws, 'numdraws', '.txt', sep=''), header=T, sep='\t')
    colnames(sample_matrix)<-indexer
} else {
    print('no sample matrix found')
    print('begin independent sampling')
    print(date())
    source('precompute_draws_unibi.r')
    print('independent sampling complete')
    print(date())
}

print('importing search methods and creating search objects')

# set the desired temperature for the system
state_temp<-1.5 

# choose temperature value that is closest to the desired
temp_ind<-which.min(abs(temp_points-state_temp))

# set initial and target states
init<-c(lambda_points[1], temp_points[temp_ind])
target<-c(lambda_points[lambda_numpoints], temp_points[temp_ind])

# set init and target lookups
init_key<-paste(init, collapse='_')
init_ind<-which(indexer==init_key)
target_key<-paste(target, collapse='_')
target_ind<-which(indexer==target_key)

# import movement definitions
source('mvmt_import.r')

# import distance evaluation methods
measure<-'asympvar'
source('dist_import_unibi.r')

# import search algorithm
source('search_import_unibi.r')

# create adjacency lists
adj_list<-vector('list', length(indexer))
names(adj_list)<-indexer
for (i in 1:length(point_list)) {
    move_bool<-apply(move_universe, 1, is_valid_move, point_list[[i]])
    valid_moves<-move_universe[move_bool,]
    neighbors<-t(apply(valid_moves, 1, move_result, point_list[[i]]))
    adj_list[[i]]<-apply(neighbors, 1, paste, collapse='_')
}

# select maximum path length of interest
# set this as the length of a path taking one-jump steps up the temperature ladder, across lambda, then back down temperature (the box solution)
max_num_nodes<-lambda_numpoints + 2*(temp_numpoints-temp_ind)

# initialize unscaled distance array
dist_matrix<-array(NA, dim=c(length(indexer), length(indexer)))

print('initialization complete')

#############
# PART 2: run search
#############
print('begin path search')
print(date())

search_out<-k_node_search(dist_matrix, max_num_nodes)
dist_matrix<-search_out[[1]]
unscaled_cost<-search_out[[2]]

print('path search complete')
print(date())

# save search result to an RData file
if (!file.exists('rdatfiles')) { dir.create('rdatfiles') }
save.image(file=paste('rdatfiles/', filID, '.RData', sep=''))

#############
# PART 3: visualize results
#############
print('visualizing results')
print(date())
source('k_node_viz_unibi.r')

print('R script complete')
print(date())
