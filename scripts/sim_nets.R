#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=T)

N = NULL
out_path = "."
prefix   =  "sim_"
ncores   =  1
max_iter = 100


help_msg <- function(){
cat(
"
Simulate level-1 phylogenetic networks

Usage: sim_networks.R N 
        --out_path out_path
        --prefix prefix
        --ncores ncores
        --max_iter max_iter

Required arguments:
  N: int; number of tips on the simulated networks
  
Optional arguments:
  --out_path: string; path to store simulated networks. (default: .)
  --prefix: string; prefix for files  (default: sim_)
  --ncores: int; number of cores (default: 1)
  --max_iter: int; this number is multiplied by ~4 so that from all 
            simulated networks we make sure that we recover max_iter number of
            level-1 networks (default: 100)
")
quit(status = 1)
}


if( length(args) < 5){
  help_msg()
}


for(i in 1:length(args)){
  
  if( i == 1 && !startsWith(args[i], '--') ){
    N = as.integer(args[i])
    
  }else if( args[i] == '--out_path' ){
    out_path = args[i + 1]
    
  }else if( args[i] == '--prefix' ){
    prefix = args[i + 1]
    
  }else if( args[i] == '--ncores' ){
    ncores   =  as.integer(args[i + 1])
    
  }else if( args[i] == '--max_iter' ){
    max_iter = as.integer(args[i + 1])
    
  }else if( args[i] == '--help' || args[i] == '-h'){
    help_msg()
  }
}

if (is.null(N) == T){
  help_msg()
}


# importing libraries
library(SiPhyNetwork)
library(parallel)


inheritance.fxn <- make.beta.draw(10,10)


#We also want to set the proportion of each type of hybrid event
hybrid_proportions <-c(0.99,  ##Lineage Generative
                       0.15, ##Lineage Degenerative
                       0.15) ##Lineage Neutral




# N = 200
nu = exp(1.6)/(N*(N-1)/2) # level 1 nets

net_names <- function(net){
  # net$tip.label <- as.numeric(sub('t', '', net$tip.label))
  net_string <- write.net( net, '')
  # gsub('t([0-9]+)','\\1', net_string)
  all_splits = strsplit( net_string, split = "(\\(|\\)|\\:|\\,|\\;)" )


  taxa_shuffle = sample(1:N, N, replace = F)
  spps = 1
  for(i in all_splits[[1]]){
    if(grepl('^t', i)){
      net_string =  sub(i, taxa_shuffle[spps], net_string)
      spps = spps + 1
    }
  }
  return(net_string)
}


i = 1
# while (T) {
start.time <- Sys.time()
suppressWarnings({
mclapply(i:as.integer(max_iter*4.6), function(i){

  all_files = Sys.glob( file.path(out_path, paste0(prefix, "*.txt") ) )

  if(length(all_files) > max_iter){
    # ('max iterations achieved')
    return(NULL)
  }

  ssa_nets<-sim.bdh.taxa.ssa(n=N,
                             numbsim=1,
                             lambda=1,
                             mu=0.2,
                             nu=nu,
                             hybprops = hybrid_proportions,
                             hyb.inher.fxn = inheritance.fxn,
                             complete = F) # ultrametric nets

  net = ssa_nets[[1]]

  # plot(net)
  if (  is.list(net) ){
    net_level = SiPhyNetwork::getNetworkLevel(net)

    # plot(net)

    if(net_level == 1){
      out_name = paste0(out_path, "/", prefix, i, '.txt')
      net_str = net_names(net)

      cat(net_str, file = out_name)
      return(net_level)
    }

  }},
  mc.cores = ncores, mc.cleanup = T) -> res
})


end.time <- Sys.time()
time.taken <- end.time - start.time


n_level_1 <- length(unlist(res))
cat(paste(n_level_1, 'level-1 random networks done: ', time.taken, 'seconds\n'))


