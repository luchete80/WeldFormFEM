class mesh

#Data
lx=1.
ly=1.
nex=100
ney=100
#***** MESH**********
dx=lx/nex
dy=ly/ney
numel=nex*ney
vn=zeros(8) #Dof connectivity
numnodes=(nex+1)*(ney+1)
node=zeros((numnodes, 2))
elnodes=zeros((numel, 4))#Connectivity
elnodes.astype(int)
#print (node)
#Mesh generation 
y=0.
n=0
for nx in range (ney+1):
    x=0.
    for ny in range (nex+1):
        #print("Node: "+str(x)+" , "+str(y))
        #this->node.push_back(Node(n,x,y,0.0))
        node[n]=[x,y]
        x=x+dx
        n=n+1
    y=y+dy
    
print(node)
#Connectivity
e=0
for ey in range (ney):
    for ex in range (nex):
        elnodes[e]=[(nex+1)*(ey+1)+ex+1,
                    (nex+1)*(ey+1)+ex,
                    (nex+1)*ey + ex,
                    (nex+1)*ey + ex+1]
        e=e+1
print(elnodes)