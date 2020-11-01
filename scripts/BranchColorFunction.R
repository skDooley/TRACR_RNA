
####### ML ASE Majority Rules Counting Transition Function - Divtree Input

MajorityRules <- function(phy, tipstates, init_input = 0.1) {
  
  if(anyNA(names(tipstates))){
    stop(paste("Object tipstates must be a named vector"))
  }
  if(class(tipstates) == "factor"){
    stop(paste("Tipstates vector must be of class character or integer"))
  }
  if(anyNA(match(phy$tip.label,names(tipstates)))){
    stop(paste("Some tips of the phylogeny do not have corresponding tipstates"))
  }
  if(anyNA(match(names(tipstates), phy$tip.label))){
    stop(paste("Some tips included in tipstates are not represented on the phylogeny"))
  }

n <- length(phy$tip.label)
traitlevels <- length(levels(as.factor(tipstates))) # calculating how many states possible for discrete variable

trait.m <- matrix(0, ncol = traitlevels, nrow = n) 
for (i in 1:n){trait.m[i,tipstates[i]] <- 1} # making a matrix of the tip labels to match ASE format

library(diversitree)
lik<-make.mkn(phy,tipstates,k=traitlevels) # constructing an MK model based on the data (MK10 bc 10 states)
init <- rep(init_input, traitlevels^2-traitlevels) # set up initial values to find best fit
fitted.ard<-find.mle(lik,x.init=init) # finding the max likelihood estimate of the data
  
asr.ard<-t(asr.marginal(lik,pars=fitted.ard$par)) # extracting the marginal likelihoods of each node (t() transposes the output to be in the right format)
colnames(asr.ard) <- levels(as.factor(tipstates)) # labeling the columns according to the tip states

divtreedata <- rbind(trait.m, asr.ard) # adding the trait matrix and the ancestral state estimations into one matrix
row.names(divtreedata) <- c(names(tipstates),(n+1):(n*2-1)) # labeling the nodes according to their node labels when plotted

node.numbers <- (length(phy$tip.label)+1):(length(phy$tip.label)*2-1) # numbering the nodes in the object phy to match the above line of code
phy$node.labels <- node.numbers

edge.mat <- phy$edge # this extraction is for the following 2 lines of code to manipulate in order to make sure tip labels are named appropriately (base and plotting systems don't align well)
tiplabel_placeholder <- phy$tip.label[edge.mat[which(edge.mat[,2] < (n+1)),2]] # creating a reordered vector of tiplabels so that they match the order in which they appear in the edge mat
edge.mat[which(edge.mat[,2] < (n+1)),2] <- tiplabel_placeholder  # relabeling tip numbers as tip names, as reordered in the line above

nodes.simp <- apply(as.data.frame(divtreedata), 1, function(x){ order(x)[ncol(divtreedata)]}) # simplifying each node to match the majority rule (just for coloring branches)
edge.mat.2 <- matrix(nodes.simp[match(edge.mat, names(nodes.simp))],ncol=2, byrow = F) # ordering these simplified nodes appropriately

transmat <- matrix(NA, nrow = ncol(divtreedata), ncol = ncol(divtreedata)) # making a transition matrix (counting the number of state transitions)
colnames(transmat) <- row.names(transmat) <- colnames(divtreedata)

for (i in 1:traitlevels) {
  for(j in 1:traitlevels){
    transmat[i,j] <- length(which(edge.mat.2[,1] == i & edge.mat.2[,2] == j)) # for each cell in the transmat, count the number of transitions
  }
}

Out <- list(transmat, edge.mat.2, asr.ard)
names(Out) <- c('Transition_Counts', 'Node_States_for_Each_Branch', "Anc_States")
Out
}