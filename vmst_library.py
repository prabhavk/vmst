# Copyright 2014-2025 Prabhav Kalaghatgi

#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
# 
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.
    
import random as r
import copy 
import math as m
import igraph as ig
import os 
import re
import subprocess as sub

########################  CLASS DECLARATIONS   ##############################
class Vertex:
    def __init__(self,name):
    # union-find
        self.rep = self
        self.comp = name
    # union-find and MLVMST 
        self.rank = 0
    # MLVMST
        self.deltaMax = 0
    # general graph
        self.name = name
        self.degree=0
        self.neighbors = []
    def AddNeighbor(self,neighbor):
        self.neighbors.append(neighbor)
        self.degree+=1
    def RemoveNeighbor(self,neighbor):
        self.neighbors.remove(neighbor)
        self.degree-=1        
    def GetNeighbors(self):
        return self.neighbors

def Find(u):
    if u.rep != u:
        u.rep = Find(u.rep)
    return u.rep

def Union(u,v):
    repOfu = Find(u)
    repOfv = Find(v)
    if repOfu == repOfv: 
        pass
    else:
        if repOfu.rank < repOfv.rank:
            repOfu.rep = repOfv
            repOfv.comp = (repOfv.comp,repOfu.comp) 
        elif repOfu.rank > repOfv.rank:
            repOfv.rep = repOfu
            repOfu.comp = (repOfv.comp,repOfu.comp)
        else:
            repOfu.rep = repOfv
            repOfv.comp = (repOfv.comp,repOfu.comp)        
            repOfv.rank+=1
                     
class Graph:
    def __init__(self):
        self.vertices={}
        self.edgeWeights={}
    def AddVertex(self,v_name):
        self.vertices[v_name] = Vertex(v_name)
    def GetVertex(self,v_name):
        return self.vertices[v_name]
    def ContainsVertex(self,v_name):
        return v_name in self.vertices
    def RemoveVertex(self,v_name):
        if v_name in self.vertices:
            v = self.vertices[v_name]
            for u in v.neighbors:
                del self.edgeWeights[tuple(sorted((u.name,v_name)))]
                u.RemoveNeighbor(v)
            del self.vertices[v_name]
        else:
            return 'Error: The vertex that is to be removed is not in the graph'             
    def AddEdge(self,u_name,v_name,w):
        if u_name not in self.vertices:
            self.AddVertex(u_name)
        if v_name not in self.vertices:
            self.AddVertex(v_name)
        u = self.vertices[u_name]
        v = self.vertices[v_name]
        u.AddNeighbor(v)
        v.AddNeighbor(u)
        self.edgeWeights[tuple(sorted((u_name,v_name)))] = w
    def ContainsEdge(self,u_name,v_name):
        return tuple(sorted((u_name,v_name))) in self.edgeWeights.keys()
    def GetEdgeWeight(self,u_name,v_name):
        return self.edgeWeights[tuple(sorted((u_name,v_name)))]
    def RemoveEdge(self,u_name,v_name):
        del self.edgeWeights[tuple(sorted((u_name,v_name)))]
        u = self.vertices[u_name]
        v = self.vertices[v_name]
        u.RemoveNeighbor(v)
        v.RemoveNeighbor(u)
    def GetLeaves(self):
        leaves=[]
        for v in self.vertices.values():
            if v.degree == 1:
                leaves.append(v)
        return leaves
    def GetVerticesInComponent(self,u_name):
        u = self.vertices[u_name]
        repOfU = Find(u)
        return repOfU.comp
    def WriteToFile(self,fileName):
        edgeListFile = open(fileName,'w')
        for u_name, v_name in self.edgeWeights.keys():
            edgeListFile.write(u_name+'\t'+v_name+'\t'+str(self.edgeWeights[tuple(sorted((u_name,v_name)))])+'\n')
        edgeListFile.close()

########################  FILE I/O ############################

def WriteDistances(distances,distancesFileName,orderedVertices=''):
    if orderedVertices=='':
        vertexNameList = [x[0] for x in distances.keys()]
        vertexNameList+= [x[1] for x in distances.keys()]
        vertexNameList = list(set(vertexNameList))
        vertexNameList.sort()
    else: vertexNameList = orderedVertices
    distanceFile = open(distancesFileName,'w')
    for i in range(len(vertexNameList)):
        for j in range(i+1,len(vertexNameList)):
            distanceFile.write(vertexNameList[i]+'\t'+vertexNameList[j]+'\t'+str(distances[tuple(sorted([vertexNameList[i],vertexNameList[j]]))])+'\n')
    distanceFile.close()

def ReadTree(treeFileName,treeFormat='edgeList',experimentName='test'):
    if treeFormat =='edgeList':
        edgeListFile = open(treeFileName,'r')
        T = ig.Graph(0)
        T.vs["name"]=[]
        for line in edgeListFile:
            if len(line.strip().split('\t'))==2:
                length=0
                vertex_i,vertex_j = line.strip().split('\t')
            else:
                vertex_i,vertex_j,length = line.strip().split('\t')
                
            for vertex in [vertex_i,vertex_j]:
                if vertex not in T.vs["name"]:
                    T.add_vertices(1)
                    T.vs[T.vcount()-1]["name"]=vertex
            vertexIndex_i = (T.vs["name"].index(vertex_i))
            vertexIndex_j = (T.vs["name"].index(vertex_j))
            T.add_edges([(vertexIndex_i,vertexIndex_j)])
            T.es[T.ecount()-1]["length"] = float(length)
        edgeListFile.close()
    elif treeFormat =='newick':
        devnull=open(os.devnull,'w')
        if os.path.isdir('/TL/euresist_phylodynamics/work/Projects/mstBasedPhylogenetics/scripts/'):
            scriptPath = '/TL/euresist_phylodynamics/work/Projects/mstBasedPhylogenetics/scripts/'
        elif os.path.isdir('/local/home/pk/Projects/mstBasedPhylogenetics/scripts/'):
            scriptPath = '/local/home/pk/Projects/mstBasedPhylogenetics/scripts/'
        
        tempTreeFileName=treeFileName+'.tempTree'
        RCommandForParsingTrees='Rscript\t'+scriptPath+'parseNewickTree.R\t'+treeFileName+'\t'+tempTreeFileName
        sub.call(RCommandForParsingTrees,stdout=devnull,shell=True)
        T = ReadTree(tempTreeFileName,'edgeList')
        sub.call('rm '+tempTreeFileName,stdout=devnull,shell=True)
        devnull.close()
    return T

def WriteTree(T,outputFileName,fileFormat='edgeList'):
    treeCopy = T.copy()
    outputFile = open(outputFileName,'w')
    if fileFormat == 'edgeList':
        for e in T.es:
            i,j = e.tuple
            outputFile.write(T.vs[i]["name"]+'\t'+T.vs[j]["name"]+'\t'+str(e["length"])+'\n')
    elif fileFormat =='newick':
        leafLabeledTree = ConvertGenerallyLabeledTreeToLeafLabeledTree(treeCopy)
        binaryTree = ConvertMultifurcatingTreeToBifurcatingTree(leafLabeledTree)
        treeInNewickFormat = GetNewickLabelOfLeafLabeledTree(binaryTree)
        treeInNewickFormat = re.sub("hiddenVertex_?l\d+-l\d+-l\d+_\d.\d+_\d*", "", treeInNewickFormat)
        treeInNewickFormat = re.sub("hiddenVertex_?l\d+-l\d+-l\d+", "", treeInNewickFormat)
        treeInNewickFormat = re.sub("pseudoHiddenVertexT?_?\d+_?\d*", "", treeInNewickFormat)
        outputFile.write(treeInNewickFormat)
    outputFile.close()

########################  TREE FUNCTIONS ############################

def GetNJTree(distanceMatrix,internalVertexIteration=1):
    fullvertList = [x[0] for x in distanceMatrix.keys()]
    fullvertList+= [x[1] for x in distanceMatrix.keys()]
    vertexNameList = list(set(fullvertList))
    vertexNameList.sort()
    T = ig.Graph()
    n = len(vertexNameList)
    T.add_vertices(n-T.vcount())
    T.vs["name"] = vertexNameList
    nodeCounter=1
    R={}
    for i in vertexNameList:
        R[i]=float(0)              
    for i in vertexNameList:
        for j in vertexNameList[vertexNameList.index(i)+1:]:
            dist = distanceMatrix[(i,j)]
            R[i]+=dist
            R[j]+=dist
    for i in vertexNameList:
        R[i]/=(n-2)
    while(len(R)>3):
        neighborDist = 10000
        for i in vertexNameList:
            for j in vertexNameList[vertexNameList.index(i)+1:]:
                neighborDist_current = distanceMatrix[(i,j)]-R[i]-R[j]
                if neighborDist_current + 10**-10 < neighborDist:
                    neighborDist = neighborDist_current
                    i_selected = i
                    j_selected = j
        newNode = 'hiddenVertex_'+str(internalVertexIteration)+'_'+str(nodeCounter)
        #print(newNode)
        T.add_vertices(1)
        T.vs[T.vcount()-1]["name"] = newNode
        nodeCounter+=1
        vertexNameList.remove(i_selected)
        vertexNameList.remove(j_selected)
        index = GetInsertIndex(vertexNameList, n-2, newNode)
        n-=1
        vertexNameList.insert(index, newNode)
        len_i = 0.5*(distanceMatrix[(i_selected,j_selected)]+R[i_selected]-R[j_selected])
        len_j = 0.5*(distanceMatrix[(i_selected,j_selected)]+R[j_selected]-R[i_selected])
        if len_i < 0: 
            len_j-=len_i
            len_i = 0
        if len_j < 0: 
            len_i-=len_j
            len_j = 0
        T.add_edges([(T.vs["name"].index(newNode),T.vs["name"].index(i_selected))])
        T.es[T.get_eid(T.vs["name"].index(newNode), T.vs["name"].index(i_selected))]["length"] = len_i
        T.add_edges([(T.vs["name"].index(newNode),T.vs["name"].index(j_selected))])
        T.es[T.get_eid(T.vs["name"].index(newNode), T.vs["name"].index(j_selected))]["length"] = len_j
        del R[i_selected]
        del R[j_selected]
        R[newNode]=0
        for vert in vertexNameList[:index]+vertexNameList[index+1:]:
            if vert < i_selected:
                newDist = distanceMatrix[(vert,i_selected)] + distanceMatrix[(vert,j_selected)]
                newDist -= distanceMatrix[(i_selected,j_selected)]
                newDist*=0.5
                R[vert] = (R[vert]*(n-1)-distanceMatrix[(vert,i_selected)]-distanceMatrix[(vert,j_selected)] + newDist)/(n-2)
                del distanceMatrix[(vert,i_selected)]
                del distanceMatrix[(vert,j_selected)]    
            elif j_selected < vert:
                newDist = distanceMatrix[(i_selected,vert)] + distanceMatrix[(j_selected,vert)]
                newDist -= distanceMatrix[(i_selected,j_selected)]
                newDist*=0.5
                R[vert] = (R[vert]*(n-1)-distanceMatrix[(i_selected,vert)]-distanceMatrix[(j_selected,vert)] + newDist)/(n-2)
                del distanceMatrix[(i_selected,vert)]
                del distanceMatrix[(j_selected,vert)]
            else:
                newDist = distanceMatrix[(i_selected,vert)] + distanceMatrix[(vert,j_selected)]
                newDist -= distanceMatrix[(i_selected,j_selected)]
                newDist*=0.5
                R[vert] = (R[vert]*(n-1)-distanceMatrix[(i_selected,vert)]-distanceMatrix[(vert,j_selected)] + newDist)/(n-2)
                del distanceMatrix[(i_selected,vert)]
                del distanceMatrix[(vert,j_selected)]
            if vert < newNode:
                distanceMatrix[(vert,newNode)] = newDist
            else:
                distanceMatrix[(newNode,vert)] = newDist
            R[newNode]+=newDist
        R[newNode]/=(n-2)
        del distanceMatrix[(i_selected,j_selected)]    
    v0,v1,v2 = vertexNameList
    d_01 = distanceMatrix[v0,v1] 
    d_02 = distanceMatrix[v0,v2]
    d_12 = distanceMatrix[v1,v2]
    newNode = 'hiddenVertex_'+str(internalVertexIteration)+'_'+str(nodeCounter)
    T.add_vertices(1)
    T.vs[T.vcount()-1]["name"] = newNode
    d_v0_h = 0.5*(d_01+d_02-d_12)
    d_v1_h = 0.5*(d_01+d_12-d_02)
    d_v2_h = 0.5*(d_02+d_12-d_01)
    if min(d_v0_h,d_v1_h,d_v2_h) < 0:
        if d_v0_h < 0:
            d_v1_h-=d_v0_h/2
            d_v2_h-=d_v0_h/2
            d_v0_h=0
        elif d_v1_h <0:
            d_v0_h-=d_v1_h/2
            d_v2_h-=d_v1_h/2
            d_v1_h=0
        else:
            d_v0_h-=d_v2_h/2
            d_v1_h-=d_v2_h/2
            d_v2_h=0
    T.add_edges([(T.vs["name"].index(v0),T.vs["name"].index(newNode))])
    T.es[T.get_eid(T.vs["name"].index(v0),T.vs["name"].index(newNode))]["length"]=d_v0_h
    T.add_edges([(T.vs["name"].index(v1),T.vs["name"].index(newNode))])
    T.es[T.get_eid(T.vs["name"].index(v1),T.vs["name"].index(newNode))]["length"]=d_v1_h
    T.add_edges([(T.vs["name"].index(v2),T.vs["name"].index(newNode))])
    T.es[T.get_eid(T.vs["name"].index(v2),T.vs["name"].index(newNode))]["length"]=d_v2_h
    return(T)
    


def ContractLatentIncidentShortBranches(T,threshold=10**-5):
    removedVertex = True
    while removedVertex==True:
        if min(T.es["length"])>threshold:
            return (T)
        removedVertex=False
        edgeIndices = range(0,T.ecount())
        r.shuffle(edgeIndices)
        for i in edgeIndices:
            vert_i, vert_j = T.es[i].tuple
            if T.es[i]["length"]< threshold:
                if T.vs[vert_i]["name"].startswith('h'):
                    neighbors_i = ig.Graph.neighbors(T,vert_i)
                    for neighbor in list(set(neighbors_i)-set([vert_j])):
                        T.add_edges([(vert_j,neighbor)])
                        lengthToVert_i = T.es[T.get_eid(vert_i,neighbor)]["length"]
                        T.es[T.get_eid(vert_j,neighbor)]["length"] = lengthToVert_i 
                        T.delete_edges([(vert_i,neighbor)])
                    T.delete_vertices(vert_i)
                    removedVertex=True
                    break
                elif T.vs[vert_j]["name"].startswith('h'):
                    neighbors_j = ig.Graph.neighbors(T,vert_j)
                    for neighbor in list(set(neighbors_j)-set([vert_i])):
                        T.add_edges([(vert_i,neighbor)])
                        lengthToVert_j = T.es[T.get_eid(vert_j,neighbor)]["length"]
                        T.es[T.get_eid(vert_i,neighbor)]["length"] = lengthToVert_j 
                        T.delete_edges([(vert_j,neighbor)])
                    T.delete_vertices(vert_j)
                    removedVertex=True
                    break
    return (T)


def GetDistanceMatrixFromTree(vertexNameList,T,idType='name'):
    # keys are index in vertexNameList
    distanceMatrix={}
    if idType =='index':
        for i in vertexNameList:
            pathsFromi = T.get_all_shortest_paths(i)
            for j in vertexNameList[vertexNameList.index(i)+1:]:
                path_ij = pathsFromi[j]
                pathLength=0
                for pathPos in range(0,len(path_ij)-1):
                    pathLength+=T.es[T.get_eid(path_ij[pathPos],path_ij[pathPos+1])]["length"]
                if i<j:
                    distanceMatrix[i,j] = pathLength
                else:
                    distanceMatrix[j,i] = pathLength
    elif idType=='name':   
        for i in vertexNameList:
            pathsFromi = T.get_all_shortest_paths(T.vs["name"].index(i))
            for j in vertexNameList[vertexNameList.index(i)+1:]:
                path_ij = pathsFromi[T.vs["name"].index(j)]
                pathLength=0
                for pathPos in range(0,len(path_ij)-1):
                    pathLength+=T.es[T.get_eid(path_ij[pathPos],path_ij[pathPos+1])]["length"]
                if i<j:
                    distanceMatrix[i,j] = pathLength
                else:
                    distanceMatrix[j,i] = pathLength
    return(distanceMatrix)

def GetRFDist(treeTrueOrig,treeEstimateOrig):
    treeTrue = copy.deepcopy(treeTrueOrig)
    treeEstimate = copy.deepcopy(treeEstimateOrig)
    for v in range(0,treeTrue.vcount()):
        if treeTrue.vs[v]["name"].startswith('hiddenVertex'):
            treeTrue.vs[v]["name"]="hiddenVertex"
    for v in range(0,treeEstimate.vcount()):
        if treeEstimate.vs[v]["name"].startswith('hiddenVertex'):
            treeEstimate.vs[v]["name"]="hiddenVertex"
    splitListTrue=[]
    for e in treeTrue.es:
        graphCopy = copy.copy(treeTrue)
        graphCopy.delete_edges([(e.tuple[0],e.tuple[1])])
        graphLets = graphCopy.decompose()
        vertList1 = graphLets[0].vs["name"]
        vertList2 = graphLets[1].vs["name"]
        removeL = 'hiddenVertex' in vertList1
        while removeL:
            vertList1.remove('hiddenVertex')
            removeL = 'hiddenVertex' in vertList1
        removeL = 'hiddenVertex' in vertList2
        while removeL:
            vertList2.remove('hiddenVertex')
            removeL = 'hiddenVertex' in vertList2
        vertList1.sort()
        vertList2.sort()
        if vertList1 < vertList2:
            splitListTrue.append(''.join(vertList1) + '-' + ''.join(vertList2))
        else:
            splitListTrue.append(''.join(vertList2) + '-' + ''.join(vertList1))
    splitListEstimate=[]
    for e in treeEstimate.es:
        graphCopy = copy.copy(treeEstimate)
        graphCopy.delete_edges([(e.tuple[0],e.tuple[1])])
        graphLets = graphCopy.decompose()
        vertList1 = graphLets[0].vs["name"]
        vertList2 = graphLets[1].vs["name"]
        removeL = 'hiddenVertex' in vertList1
        while removeL:
            vertList1.remove('hiddenVertex')
            removeL = 'hiddenVertex' in vertList1
        removeL = 'hiddenVertex' in vertList2
        while removeL:
            vertList2.remove('hiddenVertex')
            removeL = 'hiddenVertex' in vertList2
        vertList1.sort()
        vertList2.sort()
        if vertList1 < vertList2:
            splitListEstimate.append(''.join(vertList1) + '-' + ''.join(vertList2))
        else:
            splitListEstimate.append(''.join(vertList2) + '-' + ''.join(vertList1))
    jaccardIndex = len(set(splitListEstimate)&set(splitListTrue))/float(len((set(splitListEstimate)|set(splitListTrue))))
    return(1-jaccardIndex)

def BuildGraphFromEdges(edges):
    T = ig.Graph(0)
    T.vs["name"]=[]
    for vertex_i,vertex_j,length in edges:
        for vertex in [vertex_i,vertex_j]:
            if vertex not in T.vs["name"]:
                T.add_vertices(1)
                T.vs[T.vcount()-1]["name"]=vertex
        vertexIndex_i = (T.vs["name"].index(vertex_i))
        vertexIndex_j = (T.vs["name"].index(vertex_j))
        T.add_edges([(vertexIndex_i,vertexIndex_j)])
        T.es[T.ecount()-1]["length"] = float(length)
    return T

def ConvertToIgraphObject(MST_graphObj):
    edgeWeightsHash = MST_graphObj.edgeWeights
    edgesInMST = [[key[0],key[1],edgeWeightsHash[key]] for key in edgeWeightsHash.keys()]
    return BuildGraphFromEdges(edgesInMST)

def ConvertToGraphObject(graph):
    G = Graph()
    edgeNameAndWeightList = map(lambda e: graph.vs[e.tuple]["name"]+[e["length"]], graph.es)
    map(lambda e: G.AddEdge(e[0], e[1], e[2]),edgeNameAndWeightList)
    return G

def ConvertMultifurcatingTreeToBifurcatingTree(T):
    initialDegrees = T.degree()
    numberOfLatentVertices = 1
    verticesToResolve=[]
    for vertex in range(T.vcount()):
        if initialDegrees[vertex]>3:
            verticesToResolve.append(T.vs[vertex]["name"])
        if T.vs[vertex]["name"].startswith("hiddenVertex"):
            numberOfLatentVertices+=1
    for vertexName in verticesToResolve:
        while T.degree()[T.vs["name"].index(vertexName)] > 3:
            v1, v2 = T.vs[T.neighbors(T.vs["name"].index(vertexName))[0:2]]["name"]
            newNode = 'hiddenVertexT' + str(numberOfLatentVertices)
            T.add_vertices(1)
            T.vs[T.vcount()-1]["name"]=newNode
            d_v0_h = T.es[T.get_eid(T.vs["name"].index(v1),T.vs["name"].index(vertexName))]["length"] 
            T.add_edges([(T.vs["name"].index(v1),T.vs["name"].index(newNode))])
            T.es[T.get_eid(T.vs["name"].index(v1),T.vs["name"].index(newNode))]["length"] = d_v0_h
            d_v1_h = T.es[T.get_eid(T.vs["name"].index(v2),T.vs["name"].index(vertexName))]["length"]
            T.add_edges([(T.vs["name"].index(v2),T.vs["name"].index(newNode))])
            T.es[T.get_eid(T.vs["name"].index(v2),T.vs["name"].index(newNode))]["length"] = d_v1_h
            T.add_edges([(T.vs["name"].index(vertexName),T.vs["name"].index(newNode))])
            T.es[T.get_eid(T.vs["name"].index(newNode),T.vs["name"].index(vertexName))]["length"] = 0
            T.delete_edges([(T.vs["name"].index(v1),T.vs["name"].index(vertexName))])
            T.delete_edges([(T.vs["name"].index(v2),T.vs["name"].index(vertexName))])
            numberOfLatentVertices+=1
    return T
 
def ConvertGenerallyLabeledTreeToLeafLabeledTree(T):
    degrees = T.degree()
    internalLabeledVertices=[]
    largestIdOfLatentVertex=1
    for vertex in range(T.vcount()):
        if T.vs[vertex]["name"].startswith("hiddenVertex"):
            largestIdOfLatentVertex+=1
        elif degrees[vertex]>1:
            internalLabeledVertices.append(T.vs[vertex]["name"])
        else:
            pass
    for vertexName in internalLabeledVertices:
        newNode = 'hiddenVertexT'+str(largestIdOfLatentVertex)
        largestIdOfLatentVertex+=1
        T.add_vertices(1)
        T.vs[T.vcount()-1]["name"]=newNode
        neighbors = T.vs[T.neighbors(T.vs["name"].index(vertexName))]["name"]
        for v in neighbors:
            d_v0_h = T.es[T.get_eid(T.vs["name"].index(v),T.vs["name"].index(vertexName))]["length"] 
            T.add_edges([(T.vs["name"].index(v),T.vs["name"].index(newNode))])
            T.es[T.get_eid(T.vs["name"].index(v),T.vs["name"].index(newNode))]["length"] = d_v0_h
            T.delete_edges([(T.vs["name"].index(v),T.vs["name"].index(vertexName))])
        T.add_edges([(T.vs["name"].index(vertexName),T.vs["name"].index(newNode))])
        T.es[T.get_eid(T.vs["name"].index(newNode),T.vs["name"].index(vertexName))]["length"] = 0
    return T

def GetNewickLabelOfLeafLabeledTree(T): 
    # Add root at the midpoint between a leaf and its parent
    degrees = T.degree()
    T.add_vertices(1)
    T.vs[T.vcount()-1]["name"] = 'hiddenVertex_root'
    leaf = degrees.index(1)
    parentOfLeaf = T.neighbors(leaf)[0]
    length = T.es[ig.Graph.get_eid(T,leaf,parentOfLeaf)]["length"]
    T.add_edges([(leaf,T.vcount()-1)])
    T.es[T.get_eid(leaf,T.vcount()-1)]["length"]=length/2  
    T.add_edges([(parentOfLeaf,T.vcount()-1)])
    T.es[T.get_eid(parentOfLeaf,T.vcount()-1)]["length"]=length/2
    T.delete_edges([(leaf,parentOfLeaf)])
    degrees.append(2)
    newickLabels = {} 
    orderedListOfVertices = []
    numberOfTimesVertexIsVisited = {}
    root = T.vcount()-1
    for vertex in range(0,T.vcount()):
        if  degrees[vertex]>1:
            numberOfTimesVertexIsVisited[vertex]=0
            newickLabels[vertex]='('
        else:
            orderedListOfVertices.append(vertex)
            parentOfVertex = T.get_shortest_paths(root)[vertex][-2]
            newickLabels[vertex] = T.vs[vertex]["name"]+':'+str(T.es[T.get_eid(parentOfVertex,vertex)]["length"])
    while len(orderedListOfVertices)>0:
        vertex = orderedListOfVertices[0]
        if degrees[vertex]==1:
            parentOfVertex = T.neighbors(vertex)[0]
            newickLabels[parentOfVertex]+=newickLabels[vertex]
            if numberOfTimesVertexIsVisited[parentOfVertex] == degrees[parentOfVertex]-2 and T.vs[parentOfVertex]["name"]!='hiddenVertex_root':
                newickLabels[parentOfVertex]+=')'
                orderedListOfVertices.append(parentOfVertex)
            else:
                newickLabels[parentOfVertex]+=','
                numberOfTimesVertexIsVisited[parentOfVertex]+=1
        else:
            pathFromRootToVertex = T.get_shortest_paths(root)[vertex]
            parentOfVertex = pathFromRootToVertex[-2]
            newickLabels[vertex]+= T.vs[vertex]["name"]+':'+str(T.es[T.get_eid(parentOfVertex,vertex)]["length"])
            newickLabels[parentOfVertex]+=newickLabels[vertex]
            if T.vs[parentOfVertex]["name"]=='hiddenVertex_root' and numberOfTimesVertexIsVisited[parentOfVertex] == degrees[parentOfVertex]-1:
                newickLabels[parentOfVertex]+=');'
            elif numberOfTimesVertexIsVisited[parentOfVertex] == degrees[parentOfVertex]-2:
                newickLabels[parentOfVertex]+=')'
                orderedListOfVertices.append(parentOfVertex)
            else:
                newickLabels[parentOfVertex]+=','
                numberOfTimesVertexIsVisited[parentOfVertex]+=1
        orderedListOfVertices.remove(vertex)
    return newickLabels[root]


def GetInsertIndex(sortedList,lengthOfSortedList,itemToInsert,):
        # upper and lower boundaries needed to be remembered
        i=int(m.floor(lengthOfSortedList/float(2)))
        lower = 0; upper = lengthOfSortedList-1
        while True:
                if sortedList[i] > itemToInsert:
                        if i==0:
                            return(i)
                        elif sortedList[i-1]>itemToInsert:
                            upper = i-1
                            i=int(m.floor((lower+i)/float(2)))
                        else:                            
                            return(i)
                else:
                        if i==lengthOfSortedList-1:
                            return(i+1)
                        elif sortedList[i+1] <= itemToInsert:
                            lower = i+1
                            i=int(m.ceil((upper+i)/float(2)))
                        else:
                            return(i+1)



########################  ALGORITHMS FOR IMPLEMENTING VMST, MLVMST and CLG ############################

# Input: (edge weights,ordered vertices)
# Output: edges arranged in increasing order of edge weight and vertex order
def OrderEdgesOnBasisOfEdgeWeightAndVertexOrder(edgeWeights,orderedVertices,signifDigits=8):
    listToSort=[]
    for u_name,v_name in edgeWeights.keys():
        u_rank = orderedVertices.index(u_name)
        v_rank = orderedVertices.index(v_name)
        listToSort.append([round(edgeWeights[tuple(sorted([u_name,v_name]))],8),min(u_rank,v_rank),max(u_rank,v_rank),u_name,v_name])
    listToSort.sort(key=lambda x:(x[0],x[1],x[2]))    
    sortedEdges = [(x[3],x[4],x[0]) for x in listToSort]
    return sortedEdges 

# Input edge weights
# Output: An MST
def ConstructMSTWithShuffledEdgeWeights(edgeWeights,signifDigits=8):
    edgeTuple=[]
    for u_name, v_name in edgeWeights.keys():
        w = round(edgeWeights[tuple(sorted([u_name,v_name]))],8)
        edgeTuple.append([u_name,v_name,w])
    r.shuffle(edgeTuple)
    edgeTuple.sort(key=lambda x:x[2])
    M = Graph()
    for u_name, v_name, w in edgeTuple:
        if not M.ContainsVertex(v_name):
            M.AddVertex(v_name)
        if not M.ContainsVertex(u_name):
            M.AddVertex(u_name)        
        u = M.GetVertex(u_name)
        v = M.GetVertex(v_name)
        if Find(u)!=Find(v):
            M.AddEdge(u_name, v_name, w)
            Union(u, v)
    return M

# Input (edge weights,ordered vertices)
# Output: A vertex order based MST
def ConstructVMST(edgeWeights,orderedVertices):
    sortedEdges = OrderEdgesOnBasisOfEdgeWeightAndVertexOrder(edgeWeights, orderedVertices)
#     WriteDistances(sortedEdges, projectPath+'orderedDistances', orderedVertices)
    M = Graph()
    for u_name, v_name, w in sortedEdges:
        if not M.ContainsVertex(v_name):
            M.AddVertex(v_name)
        if not M.ContainsVertex(u_name):
            M.AddVertex(u_name)        
        u = M.GetVertex(u_name)
        v = M.GetVertex(v_name)
        if Find(u)!=Find(v):
            M.AddEdge(u_name, v_name, w)
            Union(u, v)
    return M

# Input: G=(V,E,w)
# Output: The common laminar family Fc, and the MST union graph Gu
def ConstructCommonLaminarFamilyAndUnionGraph(edgeWeights,signifDigits=8):
    Fc=set()
    Gu=Graph()
    edgeTuple=[]
    for u_name, v_name in edgeWeights.keys():
        w = round(edgeWeights[tuple(sorted([u_name,v_name]))],8)
        edgeTuple.append([u_name,v_name,w])
    edgeTuple.sort(key=lambda x:x[2])
    Ew=[]
    Vw=[]
    w_prev = edgeTuple[0][2]  
    for u_name, v_name, w_curr in edgeTuple:
        Fc.update([u_name,v_name])
        if not Gu.ContainsVertex(u_name):
            Gu.AddVertex(u_name)
        if not Gu.ContainsVertex(v_name):
            Gu.AddVertex(v_name)
        u = Gu.GetVertex(u_name)
        v = Gu.GetVertex(v_name)    
        if w_curr > w_prev:
            for u,v in Ew: 
                if Find(u)!=Find(v):
                    Union(u, v) 
            Ew=[]
            for u in Vw:
                repOfu = Find(u)                
                Fc.add(repOfu.comp)                
            Vw=[]
        elif Find(u)!=Find(v):
            Gu.AddEdge(u.name, v.name, w_curr)
            Ew.append((u,v))
            Vw.append(u)
        w_prev = w_curr
    return Fc, Gu
      
# Input: G=(V,E,w)
# Output: A VMST that has the minimum no. of leaves
def ConstructMinLeavesVMST(edgeWeights):
    Fc, Gu = ConstructCommonLaminarFamilyAndUnionGraph(edgeWeights)
    Fc = sorted(Fc,reverse=True)
    deltaMax={}
    for v in Gu.vertices.values():
        deltaMax[v]=0
        N = set(u.name for u in set(v.neighbors))
        for C in Fc:
                if v.name not in C and len(N&set(C))>0:
                    deltaMax[v]=deltaMax[v]+1
                    N=N-set(C)
                
    deltaMax_list = [(v.name,deltaMax[v]) for v in deltaMax.keys()]
    deltaMax_list.sort(key=lambda x:x[1])
    orderedVertices = [x[0] for x in deltaMax_list]
    M = ConstructVMST(edgeWeights, orderedVertices)
    return M

def GetSubsetOfDistances(distances,vertexNameList):
    subsetDistances={}
    vertexNameList.sort()
    for u_name in vertexNameList:
        for v_name in vertexNameList[vertexNameList.index(u_name)+1:]:
            subsetDistances[(u_name,v_name)] = distances[(u_name,v_name)]
    return subsetDistances

def tsort(x,y):
    return tuple(sorted([x,y]))

# Input: MST, distances d
# Output: Phylogenetic tree
def CLG(M,d):
    Vm = set(M.vertices.values())
    Vm_names = set(M.vertices.keys())
    I_names = [u.name for u in list(Vm-set(M.GetLeaves()))]
    T = Graph()
    for (u_name, v_name) in M.edgeWeights.keys():
        w = M.GetEdgeWeight(u_name, v_name)
        T.AddEdge(u_name, v_name, w)
    s = {}
    for u_name in I_names:
        u = T.GetVertex(u_name)
        Vu = set(u.neighbors+[u])
        Vu_names = set([x.name for x in Vu])
        if len(Vu_names-Vm_names)>0:
            hInVu_names = list(Vu_names-Vm_names)
            for h_name in hInVu_names:
                T.RemoveEdge(u_name, h_name)
                for v_name in list(Vu_names&Vm_names):
                    d[tsort(h_name,v_name)] = d[tsort(s[h_name],v_name)]-d[tsort(h_name,s[h_name])]
            for h_1_name in hInVu_names:
                for h_2_name in hInVu_names[hInVu_names.index(h_1_name)+1:]:
                    d[tsort(h_1_name,h_2_name)] = d[tsort(s[h_1_name],s[h_2_name])] - d[tsort(h_1_name,s[h_1_name])]-d[tsort(s[h_2_name],h_2_name)]                    
        Tu=GetNJTree(GetSubsetOfDistances(d, list(Vu_names)),I_names.index(u_name)+1)
        dTu = GetDistanceMatrixFromTree(Tu.vs["name"],Tu)
        Tu = ConvertToGraphObject(Tu)
        VTu_names = set(Tu.vertices.keys())
        for h_name in VTu_names-Vm_names:
            s[h_name] = u_name
            d[tsort(h_name,u_name)] = dTu[tsort(h_name,u_name)]
        Nmu_names = [x.name for x in M.GetVertex(u_name).neighbors]
        for v_name in Nmu_names:
            if T.ContainsEdge(u_name, v_name):
                T.RemoveEdge(u_name, v_name)
        for i_name, j_name in Tu.edgeWeights.keys():
            w = Tu.GetEdgeWeight(i_name, j_name)
            T.AddEdge(i_name, j_name, w)
    return T


if False:
    # Change project path as appropriate.
    projectPath = '/home/pk/Projects/vertex-rankedMSTs/results/primatePhylogeny/'
    T = ReadTree(projectPath+'Primates_genus_zeroLengthEdgesContracted.nwk', 'newick')
    vertexNameList=[]
    for v_name in T.vs["name"]:
        if not v_name.startswith('hidden'):
            vertexNameList.append(v_name)
    
    
    orderedVertices = vertexNameList
    d = GetDistanceMatrixFromTree(orderedVertices, T)
    WriteDistances(d, projectPath+'distances')
    r.shuffle(orderedVertices)
    M = ConstructVMST(d, orderedVertices)
    print ('randomly shuffled order', GetRFDist(T, ConvertToIgraphObject(M)))
    CLGTree = ConvertToIgraphObject(CLG(M, d))
    CLGTree = ContractLatentIncidentShortBranches(CLGTree,10**-4)
    print (GetRFDist(T, CLGTree))
    
    orderedVertices.sort(reverse=False)
    d = GetDistanceMatrixFromTree(orderedVertices, T)  
    M = ConstructVMST(d, orderedVertices)
    print ('ascending order', GetRFDist(T, ConvertToIgraphObject(M)))
    CLGTree = ConvertToIgraphObject(CLG(M, d))
    CLGTree = ContractLatentIncidentShortBranches(CLGTree,10**-4)
    print (GetRFDist(T, CLGTree))
    
    orderedVertices.sort(reverse=True)
    d = GetDistanceMatrixFromTree(orderedVertices, T)
    M = ConstructVMST(d, orderedVertices)
    print ('descending order', GetRFDist(T, ConvertToIgraphObject(M)))
    CLGTree = ConvertToIgraphObject(CLG(M, d))
    CLGTree = ContractLatentIncidentShortBranches(CLGTree,10**-4)
    print (GetRFDist(T, CLGTree))
    
    # MST that is not necessarily a VMST   
    maxRFDist=0
    for _ in range(100):
        d = GetDistanceMatrixFromTree(orderedVertices, T)
        M = ConstructMSTWithShuffledEdgeWeights(d)
        CLGTree = ConvertToIgraphObject(CLG(M, d))
        CLGTree = ContractLatentIncidentShortBranches(CLGTree,10**-4)
        if GetRFDist(T, CLGTree) > maxRFDist:
            maxRFDist = GetRFDist(T, CLGTree)
            print (maxRFDist)
            CLGTree_maxRFDist = CLGTree    
    WriteTree(CLGTree_maxRFDist, projectPath+'CLGTree_maxRFDist', 'newick')