# search_import.r

# algorithm from http://stackoverflow.com/questions/1690347/shortest-path-with-a-fixed-number-of-edges 

k_node_search<-function(D, k) {
    # D[i,j] gives unscaled distance going from node i to node j
    # k is the number of nodes in the longest solution you want to know about

    # initialize cumulative best distance array, where entry i,j is the best distance to node i in j steps
    part_dist<-array(Inf, dim=c(length(indexer), k))
    part_dist[init.ind,1]<-0
    rownames(part_dist)<-indexer

    # initialize array that tracks path
    prev_mat<-array(NA, dim=c(length(indexer), k))
    rownames(prev_mat)<-indexer

    # for each path length of interest, find the best path to each node
    for (i in 2:k) {
        print(paste('searching for paths of length ', i, sep=''))

        # which nodes can even connect to the source with i-1 links?
        # get nodes with a path to it with  i-2 links
        connected_nodes<-which(part_dist[,i-1]!=Inf)

        # get a list of the neighbors of those connected nodes
        query_nodes<-unique(unlist(adj_list[connected_nodes]))

        # for each of those nodes, look at each of its neighbors
        for (j in 1:length(query_nodes)) {
            j_ind<-which(indexer==query_nodes[j])
            j_neighbors<-adj_list[[j_ind]]

            # for each neighbor, get the cost to the node through that neighbor
            for (l in 1:length(j_neighbors)) {

                # get cost up to neighbor
                neigh_ind<-which(indexer==j_neighbors[l])
                neigh_part<-part_dist[neigh_ind, i-1]

                # get cost from neighbor to node
                # if we haven't calculated it before, do so
                if (is.na(D[j_ind, neigh_ind])) {
                    coords1<-as.numeric(unlist(strsplit(indexer[j_ind], split='_')))
                    coords2<-as.numeric(unlist(strsplit(j_neighbors[l], split='_')))
                    neigh_rest<-dist.move(coords1, coords2)
                    D[j_ind, neigh_ind]<-neigh_rest
                    D[neigh_ind, j_ind]<-neigh_rest
                # otherwise just look it up 
                } else { neigh_rest<-D[j_ind, neigh_ind] }

                # get total cost to node through neighbor
                neigh_tot<-neigh_part + neigh_rest

                # if this is best cost to node, update
                if (neigh_tot < part_dist[j_ind, i]) { 
                    part_dist[j_ind, i]<-neigh_tot
                    prev_mat[j_ind, i]<-j_neighbors[l]
                }
            } # close loop over neighbors
        } # close loop over nodes being queried
    } # close loop over prospective path lengths

    # return the updated distance matrix, the cost matrix and the path tracer
    return(list(D, part_dist, prev_mat))
}

path_soln<-function(trace_mat, num_node) {
    seq_out<-array(NA, dim=c(num_node, 2))
    seq_out[num_node,]<-target

    trace_index<-which(indexer==paste(target, collapse='_'))

    for (i in num_node:2) {
        next1<-trace_mat[trace_index, i]
        seq_out[i-1,]<-as.numeric(unlist(strsplit(next1, split='_')))
        trace_index<-which(indexer==next1)
    }

    rownames(seq_out)<-c(1:nrow(seq_out))
    colnames(seq_out)<-c('lambda','temp')
    seq_out<-as.data.frame(seq_out)
    return(seq_out)
}