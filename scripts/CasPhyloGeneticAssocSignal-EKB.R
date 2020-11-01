#setwd("/Users/shanedooley/Documents/Research/CRISPR_CAS/data")
setwd("~/Documents/School/Projects/2020-Dooley_et_al/")


## 1: Read data, tree, and prune/match one to the other
library(corHMM)
library(geiger) #also loads ape
library(phytools)
library(diversitree)
source("DiscPhySignal.R")
source("BranchColorFunction.R")
casTree = read.tree("Cas9-Like-clustered.faa.tree",tree.names = TRUE)
structureData = read.csv('StructureNumber.csv', row.names=1, header=TRUE)
data.pruned = treedata(phy=casTree,data = structureData, warnings=FALSE)
fullPruned = data.pruned$phy
tree=fullPruned
Verified <- row.names(structureData)[which(structureData$Is.Ref == "TRUE")]


# Subsample the tree, keeping all tips that have been verified
# ALL THIS HAS BEEN COMMENTED OUT AFTER AN APPROPRIATE TREE WAS IDENTIFIED: Subtree_that_works-All_Verified_Included.tre 
      # This was done because it took many iterations of subsampling to get the code to run where we constrained 
      # our analyses to keeping all the verified tips on the tree. There are certainly other trees that fit the bill,
      # but the one in the phylogeny plot is the Subtree_that_works-All_Verified_Included.tre file mentioned above

#SubSampleTree = function(tree,perc){ # this function was changed to make sure that all Verified tips were kept
#  m<-round(((1-perc)*length(tree$tip)))
#  drops<-sample(tree$tip.label[-which(tree$tip.label %in% Verified)])[1:m] # do not drop any Verified tips
#  subTree<-drop.tip(tree,drops)
#  return(subTree)
#}
#subSamplePerc = 25/100 # 
#subtree = SubSampleTree(fullPruned,subSamplePerc)

subtree <- read.tree("Subtree_that_works-All_Verified_Included.tre")
trait = structureData[subtree$tip.label,]

# 3: The analysis for test 1
tree.ult <- force.ultrametric(subtree)
y <- trait$Struct
names(y) <- row.names(trait)

length(levels(as.factor(trait$Struct)))==10 # want TRUE
      # This verifies that the subtree selected has at least one tip in each of the 10 states
   
Results <- MajorityRules(tree.ult, y, init_input = .1) 
    # This uses the BranchColorFunction script to run the ancestral state estimations, extracting the most likely state of each node
          # which was used to color the branches below
    # Anc_States (below) determined using diversitree


# 2: Parameters and functions for the analysis
lambda0 <- 0.1   #rate parameter of the proposal 
se      <- 0.5   #standard deviation of the proposal
sim     <- 100000     #number of iterations
thin    <- 10    #we kept only each 10th iterate 
burn    <- 100   #100 iterates are burned-in

ar = Results$Anc_States # aceStats$lik.anc
x  <- nentropy(ar)
mc1    <- emcmc(rexp(1),rexp(1),x,lambda0,se,sim,thin,burn)
mc2    <- emcmc(rexp(1),rexp(1),x,lambda0,se,sim,thin,burn)
mchain <- rbind(mc1,mc2)
deltaA <- mean(mchain[,2]/mchain[,1])
print(deltaA)

# 4. DRAW THE TREE...

library(RColorBrewer)
pie_cols <- brewer.pal(10, "Set3") # selecting the color spectrum from RColorBrewer
pie_cols[2] <- "#696969" # changing the yellow color to a dark gray

edge.col <- pie_cols[Results$Node_States_for_Each_Branch[,2]] # extracting the end state of each node, with which we color each branch
label.names <- LETTERS[1:10] # naming the states

# changing positions of labels manually to avoid text overlap
shortened.ref.code.all <- matrix(unlist(strsplit(row.names(trait)[which(trait$Is.Ref =="TRUE")], "_")), ncol = 4, byrow = T)
shortened.ref.code <- shortened.ref.code.all[,3]
star_pos <- rep(" ", nrow(trait))
star_pos[which(trait$Is.Ref =="TRUE")] <- shortened.ref.code
shortened.ref.code[which(shortened.ref.code %in% c("Seq2", "Sra", "Tpu", "Cme3"))] <- rep(" ", 4)
star_pos <- rep(" ", nrow(trait))
star_pos[which(trait$Is.Ref =="TRUE")] <- shortened.ref.code

setEPS()
postscript("Subtree_colored-V1.eps")
plot(tree.ult,no.margin = TRUE,type="fan",show.tip.label=F,edge.color=edge.col,edge.width=1.8)
tiplabels(star_pos,offset=0.03,cex=.4,bg = "transparent", col="black", 
          frame="none",srt=0)
text(-.86,0.76, "Cme3", cex = .4)
text(.24, 1.113, "Tpu", cex = .4)
text(.8,.81, "Sra", cex = .4)
text(.93,.66, "Seq2", cex = .4)
legend(-1.1,1.1, label.names, fill = pie_cols, cex = .6)
dev.off()


######### 2nd Test: Adams & Collyer 2018 Evolution Phylogenetic Anova ##########
#Using partial least squares
setwd("~/Documents/School/Projects/2020-Dooley_et_al/")
library(geomorph)
library(ape)
library(geiger)
casTree = read.tree("Cas9-Like-clustered.faa.tree",tree.names = TRUE)
structureData = read.csv('StructureNumber.csv', row.names=1, header=TRUE)
data.pruned = treedata(phy=casTree,data = structureData, warnings=FALSE)
fullPruned = data.pruned$phy

trait.m <- matrix(0, ncol = 10, nrow = length(fullPruned$tip.label)) 
for (i in 1:length(fullPruned$tip.label)){trait.m[i,data.pruned$data[,1][i]] <- 1} # making a matrix of the tip labels to match ASE format

phyloCovarMatrix = vcv.phylo(fullPruned)  # converting the phylogeny into a covariance matrix
rownames(phyloCovarMatrix) <- colnames(phyloCovarMatrix) <- NULL

partLS = two.b.pls(phyloCovarMatrix,trait.m,iter=999) # two block partial least squares test
summary(partLS)


















## SCRATCH SPACE
# subtype = subtyping[tree$tip.label,]
# structureData = read.csv('tables/Metadata.csv', row.names=1, header=TRUE) #StructureNumber.csv
# library(dplyr)
# subtyping = select(structureData,-Struct)
# structureData= select(structureData,-Type)

#Add Subtype values
pie_cols2 = c('#ff0000','#00ff00','#0000ff')
msubtype <- matrix(0,ncol=3,nrow=numLeaves)
for ( i in 1:numLeaves) {
  msubtype[i,as.numeric(subtype[i])] <- 1
}
tiplabels(pie=msubtype,cex=0.10,piecol=pie_cols2)

