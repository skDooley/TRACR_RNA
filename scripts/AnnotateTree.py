from ete3 import Tree, TreeStyle, TextFace, NodeStyle, CircleFace
import sys

nstyle = {
"SCRUBBEDSPECIESCODE":"#1f73c6", "SCRUBBEDSPECIESCODE":"#8e3ce0", "SCRUBBEDSPECIESCODE":"#ff00ff", "SCRUBBEDSPECIESCODE":"#faff0a", 
"SCRUBBEDSPECIESCODE":"#f28213", "SCRUBBEDSPECIESCODE" :"#ffccff", "SCRUBBEDSPECIESCODE":"#006600", "SCRUBBEDSPECIESCODE":"#028ef9", 
"SCRUBBEDSPECIESCODE":"#cc3300", "SCRUBBEDSPECIESCODE":"#6b996b", "SCRUBBEDSPECIESCODE_short":"#00f2ff" }


referenceNodes = {
#SCRUBBEDIDS MAPPED TO SPECIES CODES
"Locus9":"L9",
"Locus14":"L14",
"Locus15":"L15",
"Locus24":"L24",
"Locus26":"L26",
"Locus27":"L27",
"Locus28":"L28",
}

def correctIDs(IDs):
	retSet = {}
	if type(IDs) == type({}):
		for id in IDs: retSet[id.replace(":","_").replace("-","")] = IDs[id]
	else:
		for id in IDs: retSet[id.replace(":","_").replace("-","")] = "#F4B266" #"#70857A"
	return retSet
def correctNWK(nwkTree):
	retDict = {}
	for leaf in nwkTree: retDict[leaf.name] = leaf.name[leaf.name.find("_")+1:]
	return retDict
class AnnotateTree:
	def __init__(self,setOfIDs,description="",tree="data/Cas9Phylogeny.nwk",colors={},highlight=set()):
#		print "Building Tree"
		CasTree = Tree(tree)
		self.treeIDs = correctNWK(CasTree)
		correctedSet = correctIDs(setOfIDs)
		numWithHit, numWithoutHit = 0, 0
		print ("Number of leaves:", len(CasTree))
		for leaf in CasTree:
#			if "Locus" not in str(leaf.name): id = "_".join(str(leaf.name).replace("'","").split("_")[:-1])
#			else: id = str(leaf.name)
			id = str(leaf.name)
			if id in referenceNodes: 
				print (leaf.name)
				leaf.add_face(CircleFace(1000,"#d3d3d3",style='circle',label={'text':referenceNodes[id],'color':'#1f73c6','fontsize':620}), 0, position="aligned")
			node_style = NodeStyle()
			#if self.treeIDs[id] not in highlight: leaf.name = ""
					
			if id in correctedSet:
				node_style["shape"] = "sphere"
				node_style["size"] = 50
				node_style["fgcolor"] = correctedSet[id]
				node_style["bgcolor"] = correctedSet[id]
				node_style["hz_line_color"] = correctedSet[id]
				node_style['hz_line_width'] = 50
				node_style['vt_line_width'] = 50
				leaf.set_style(node_style)
#				leaf.add_face(TextFace(str(leaf.name), fsize=288,fgcolor="#000000"),0)
				numWithHit+=1
			else: 
#				leaf.name = ""
				numWithoutHit += 1
			if self.treeIDs[id] in colors:
				node_style["shape"] = "sphere"
				node_style["size"] = 50
				node_style["fgcolor"] = colors[self.treeIDs[id]]
				node_style["bgcolor"] = colors[self.treeIDs[id]]
				node_style["hz_line_color"] = colors[self.treeIDs[id]]
				node_style['hz_line_width'] = 50
				node_style['vt_line_width'] = 50
				leaf.set_style(node_style)
	  
				
#		print "Number in tree", numWithHit
#		print "Missing from tree", numWithoutHit
		ts = TreeStyle()
		ts.show_leaf_name = True #len(highlight)>0
#		print "Default tree width", ts.tree_width
		ts.mode = "c"
#		ts.tree_width = 900
		ts.optimal_scale_level ="full"
#		tf = TextFace(description, fsize=288,fgcolor="#00ff00")
#		ts.title.add_face(tf, column=0)
#		ts.allow_face_overlap =True
#		print "show"
#		CasTree.show(tree_style=ts)
#		print "Rendering"
#		CasTree.render("CasTree.png",tree_style=ts) #, w=800, units="mm"
#		print "Done"
		self.tree = CasTree
		self.ts = ts
	def contains(self,id):return id in self.treeIDs
		
		


if __name__ == "__main__":
	print ("start")
	tree = AnnotateTree(set(),"","data/CasTree/nwks/Cas9-like_final.nwk")
#	tree = AnnotateTree(set(),"","SCRUBBEDID_pam_expt.nwk")
	#tree = AnnotateTree(set(),"","data/Cas9Phylogeny.nwk")
	print("generated tree")
	
	print ("Done")
'''
allCas9Tree = Tree("data/Cas9Phylogeny.nwk")

#for node in t:
#	if node.is_leaf():
#		print node
#		print dir(node)
#		break

nstyle = {SCRUBBEDSPECIESCODESMAPPEDTOCOLORS}
hits = {}
for row in DictReader(open("data/AllStructBLASTResults.txt","r"),delimiter="\t"):
	if row['Hit'] in hits:
		if float(row['LengthMatch']) > hits[row['TreeName']]: hits[row['TreeName']]=row
	else: hits[row['TreeName']]=row
for name,color in nstyle.iteritems():
	newNode = NodeStyle()
	newNode["shape"] = "sphere"
	newNode["size"] = 5500
	newNode["fgcolor"] = color
	newNode["hz_line_color"] = color
	newNode['hz_line_width'] = 5000
	
	nstyle[name] = newNode

count = 0
countMissing = 0
for leaf in allCas9Tree:
	if str(leaf.name).replace("'","") in hits: 
		leaf.set_style(nstyle[hits[str(leaf.name).replace("'","")]['Struct']])
		count+=1
#		print "Found", leaf.name
	else:
		countMissing += 1
#		leaf.add_features(color="none")
print "Number in tree", count
print "Missing from tree", countMissing

#print clade1Tree.get_ascii(attributes=["name", "color"], show_internal=False)

ts = TreeStyle()
ts.show_leaf_name = True
#ts.show_branch_length = True
ts.mode = "c"
ts.title.add_face(TextFace("Cas9 Clade 1", fsize=20), column=0)
ts.show_branch_support = False
allCas9Tree.show(tree_style=ts)
'''


'''
http://etetoolkit.org/docs/latest/tutorial/tutorial_drawing.html
http://etetoolkit.org/docs/latest/tutorial/tutorial_trees.html#reading-newick-trees
from ete3 import Tree
t =  Tree("((((((4, e), i), o),h), u), ((3, 4), (i, june)));")
# we annotate the tree using external data
colors = {"a":"red", "e":"green", "i":"yellow",
		  "o":"black", "u":"purple", "4":"green",
		  "3":"yellow", "1":"white", "5":"red",
		  "june":"yellow"}
for leaf in t:
	leaf.add_features(color=colors.get(leaf.name, "none"))
print t.get_ascii(attributes=["name", "color"], show_internal=False)

#				   /-4, green
#				/-|
#			 /-|   \-e, green
#			|  |
#		  /-|   \-i, yellow
#		 |  |
#	   /-|   \-o, black
#	  |  |
#	/-|   \-h, none
#   |  |
#   |   \-u, purple
# --|
#   |	  /-3, yellow
#   |   /-|
#   |  |   \-4, green
#	\-|
#	  |   /-i, yellow
#	   \-|
#		  \-june, yellow

'''
