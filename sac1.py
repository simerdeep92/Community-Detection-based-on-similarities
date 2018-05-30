import sys
import csv
import math
import igraph
from igraph import *
from itertools import izip
from scipy import spatial

##Return Cosine Similarity
def CosineSimilarity(v1,v2):
	return 1 - spatial.distance.cosine(v1,v2)

##find similarity matrix  - sim
def SimilartyMatrix(gra):
	sim = []	
	#d = len(headers) ## gives out the number of attributes 
	for i in gra.vs():
		temp = []
		for j in gra.vs():
			#temp.append(1.0/ DenominatorSimilarityMatrix(attrList[i],attrList[j]))
			temp.append(CosineSimilarity(i.attributes().values(),j.attributes().values()))	
		sim.append(temp)
	return sim

##Returns Attribute Similarity
def QAttr1(x,C,sim):
	result = 0
	for ver in C:
		if( x != ver):
			result = result + sim[x][ver]
	return result/(math.pow(len(C),2))

## New similarity of contracted graph using summing and normalizing the similarties of previous matrix
def NewSimilartyMatrix(newNodes,sim):
	newSim = [[0]*len(newNodes)]*len(newNodes)
	for n1,i in enumerate(newNodes):
		for n2,j in enumerate(newNodes):
			temp = 0
			count = 0
			for k in i:
				for l in j:
					temp = temp + sim[k][l]
					count = count + 1

			newSim[n1][n2] = temp/count
	return newSim	

#Phase 1 Algorithm
def phase1(membership,sim,g):
	
	maxmodularityCommunity = membership
	vertices = range (len(set(membership)))

	for loop in range(0,15):
		totalModularityChange = 0	
		for i in vertices:
			deltaCompositeModularity = 0
			maxCom = -1
			initialQNewman = g.modularity(maxmodularityCommunity) 
			#print g.modularity(maxmodularityCommunity)
			temp = membership[i]
			for j in vertices:
				
				###Check to move in different communities
				if membership[i] != membership[j] :
					#  Community of i
					Cold = [l for l, m in enumerate(membership) if m == membership[i]]
					# Move i to j's community
					membership[i] = membership[j]
					# New Community of i = j
					Cnew = [p for p, q in enumerate(membership) if q == membership[j]]
					#Calculate change in newmans modularity 				
					deltaQNewMan = g.modularity(membership) - initialQNewman
					##Return normalized change in sum of attribute similarity of Communitites
					attrSimilarity = QAttr1(i,Cnew,sim) - QAttr1(i,Cold,sim)
					finalAttrSimilarity = attrSimilarity/len(set(membership))
					## Total Modularity Change 				 
					DeltaQ = (alpha * deltaQNewMan) +  ((1-alpha) * finalAttrSimilarity)				
					## Find communitty with maximum positive gain
					if DeltaQ > deltaCompositeModularity:
						maxCom = membership[j] 
						deltaCompositeModularity = DeltaQ
					membership[i] = temp	
				else :
					deltaCompositeModularity = 0
					maxCom = -1
			## Move i to community with maximum modularity gain
			if maxCom != -1:	
				membership[i]  = maxCom				
			totalModularityChange = totalModularityChange + deltaCompositeModularity	
		##No more changes in communitities if there is no total gain in modularity
		if totalModularityChange == 0:
			break
	return membership

def main(argv):
	########Input parameter - alpha ##########
	global alpha 
	alpha = float(argv)
	global attrList 
	attrList = []
	#########Read graph using edges from file
	g = igraph.Graph()
	g = g.Read_Edgelist("data/fb_caltech_small_edgelist.txt")
	columns = defaultdict(list)
	with open('data/fb_caltech_small_attrlist.csv', 'rb') as f:
	    reader = csv.DictReader(f)	
	    for row in reader:
        	for (k, v) in row.items():
       		     columns[k].append(int(v))
	for key in columns:
    		g.vs[key] = columns[key]		
	########Assign Attribute to graph and populate attribute list
	i = 0
	f1 = open ("data/fb_caltech_small_attrlist.csv")
	reader = csv.reader(f1) 
	headers = reader.next()
	len(headers)
	for row in reader:   
		attrList.append(map(int,row))  
		for attrName,attrValue in izip(headers,row):
			g.vs[i][attrName] = float(attrValue)
		i = i+1	
	for edge in g.es():
   		 edge["weight"] = 1

	# Initialize each node to its own community
	clusters = Clustering(range(324))
	membership = clusters.membership
	
	########PHASE 1 Iteration ###########
	sim = SimilartyMatrix(g)
	phase1Membership = phase1(membership,sim,g)
	communitiesCluster = Clustering(phase1Membership)
	phase1Communities = [x for x in list(communitiesCluster) if len(x) != 0 ]
	##Writes Communities after phase 1 to file	
	with open("communities_afterphase1"+".txt", "wb") as f:
   		writer = csv.writer(f)
   		writer.writerows(phase1Communities)

	#######PHASE 2 Iteration ############
	## Simplify Membership from 0 to number of communities formed ##
	simpDict = {}
	newMemberShip = []
	count = 0 
	for i in phase1Membership:
		if i in simpDict.keys():
			newMemberShip.append(simpDict[i])
		else:
			newMemberShip.append(count);
			simpDict[i] = count
			count = count + 1
	
	## Contract the graph vertices based on phase1Communities formed ##
	g.contract_vertices(newMemberShip)
	g.simplify(multiple = True)
	## Reorder the communitites to variable newnodes
	newNodes = list(Clustering(newMemberShip))
	## Calculate new similarity matrix
	newSimilarityMatrix = NewSimilartyMatrix(newNodes,sim)
	## Run phase 1 on communities acting as a vertice 
	phase2Membership = phase1(range(len(newNodes)),newSimilarityMatrix,g) 
	##Communities after repeating phase 1 on contracted graph
	communitiesCluster = list(Clustering(phase2Membership))
	phase2Communities = [x for x in list(communitiesCluster) if len(x) != 0 ]
	## Decode new nodes for phase 2 to original vertices.
	finalCommunities = []
	for com in phase2Communities:
		temp = []
		for ver in com:
				for nodes in newNodes[ver]:
					temp.append(nodes)
		finalCommunities.append(temp)
	## Write the Coomunities data to file
	with open("communities"+".txt", "wb") as f:
   		writer = csv.writer(f)
   		writer.writerows(finalCommunities)
			
## Check for arguments and run the main algorithm
if __name__ == "__main__":
	if(len(sys.argv) == 2):  
		main(sys.argv[1])
	else :
		print "Invalid no. of arguments"
	