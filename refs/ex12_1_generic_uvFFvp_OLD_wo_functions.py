# https://docs.sympy.org/latest/tutorial/matrices.html
#Plain strain/stress 1 element example
#import math
from numpy import *
import numpy.matlib #for zeros

import array as arr

import deriv

#----------------------------
#Input Data------------------
#############################

#CASE RADIAL FLOW
# form=2
# lx=1.
# ly=20.*3.1415926/180.
# nex=2
# ney=1

form=2
lx=0.1
ly=0.01
nex=2
ney=1
plastic=0
plastic=0
numit=1
solver=2 #1:simple 2:Newton Raphson

bcnodecount=2*(ney+1)   #INLET AND CENTER 
is_bcnode_byvar=zeros((bcnodecount, 4)) #in case of disturbed flow
for n in range(ney+1):
    is_bcnode_byvar[n]=[1,1,0,0]    #Velocity and F are bc 

n=ney+1
for ni in range (ney+1):
    is_bcnode_byvar[ni+n]=[1,0,0,0] #ONLY VELOCITY IS BC

Uinit=[1.,0.,1.,0.,0.,1.] #Initial Condition, Velocity and F
#-------------
if plastic==0:
	numvars=2 #u,F
else:
	numvars=4 #u,F,F,Fvp
#BEFORE BOUNDARY CONDITIONS
#Element dimension and DOF PER VARIABLE! 
var_edof=zeros(2)
if plastic == 0:
	if form==2:
	    var_dim =[2,4]
else:
	if form==2:
		var_dim =[2,4,5,1]
	else:
		var_dim =[2,4,4,1]
var_edof=zeros(4)
ndof=0
#Formulation type and DOFs
for i in range(numvars): 
    ndof += var_dim[i]
print ("dofs: ", ndof)
for i in range(numvars):
    var_edof[i]=4*var_dim[i] 

#### BOUNDARY CONDITIONS
# BOUNDARY CONDITIONS
boundarynode=zeros(2*(ney+1))
#ROWS ARE NODES; COLS ARE ENTIRE BCs vector (all vars)
node_bc=matrix(numpy.matlib.zeros((size(boundarynode), ndof)))

dnode=(nex+1)    
i=0
for dy in range(ney+1): 
    inode=dy*dnode
    boundarynode[i]=inode
    print("i",i)
    node_bc[i]=[1.,0.,1.,0.,0.,1.]
    i+=1

inode=nex/2 #Central nodes
for dy in range(ney+1):
    boundarynode[i]=inode
    inode+=dnode
    print("i",i)
    node_bc[i]=[1.001,0.,0.,0.,0.,0.] #ONLY VELOCITY
    i+=1
    
#---------------------------------------------------
#Material properties (Table 2.1 p28)
#HSLA-65 steel with strain rate between e-3 and e-4
#---------------------------------------------------

mat_A0 =6.34e11
mat_psi=3.25
mat_m  =0.1956
mat_n  =0.06869
mat_s0 =80.0e6
mat_Q  =312.35e3
mat_R  =8.314
mat_a  =1.5
mat_h0 =3093.1e6

# print (K)
# print (X2[0,0])
# print (X2[2,1])

#--------- EXAMPLE RADIAL FLUX ----------------
ey=1.e9
nu=0.499


#***** MESH**********
dx=lx/nex
dy=ly/ney
numel=nex*ney

vnrow=zeros(20) #Dof connectivity
vncol=zeros(20) #Dof connectivity

numnodes=(nex+1)*(ney+1)
node=zeros((numnodes, 2))
vnxy=zeros((numnodes, 2))
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
 
# #Radial flow example: convert from cilindrical to cartesian
# r0=1.
# r1=2.
# for n in range(numnodes):
   # r=node[n,0]+r0
   # t=node[n,1]-ly/2.
   # node[n,0]=r*cos(t)
   # node[n,1]=r*sin(t)
   # print("Coord ",n,":",node[n,0],node[n,1])
   

# #Only for this example   
# for n in range(numnodes):
    # r=node[n,0]+r0
    # t=node[n,1]-ly/2.
    # vr=0.1*r0/r
    # vnxy[n,0]=vr*cos(t)
    # vnxy[n,1]=vr*sin(t)
# print("vxy ",n,":",vnxy[n,0],vnxy[n,1])
# #Radial flow example: convert from cilindrical to cartesian
# for n in range(numnodes):
   # r=node[n,0]+r0
   # t=node[n,1]-ly/2.
   # node[n,0]=r*cos(t)
   # node[n,1]=r*sin(t)
   # #print("Coord ",n,":",node[n,0],node[n,1])
  
    

    
print("node_bc",node_bc)
    
print("boundarynode",boundarynode)

#########################################################################################
##################################### INPUT END
 
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
#-------------------------- MESH
    
#print ("var_edof",var_edof)

edof=ndof*4
dof=ndof*numnodes
#X is combines per row to calculate J
p=1.0/1.732050807568877
gauss=[-p,p]
#Numerated as in Bathe
X2=numpy.matlib.zeros((4, 2))

#Main variables values
#Nodal values
U   =matrix(numpy.matlib.zeros((44, 1)))
UV  =matrix(numpy.matlib.zeros((8, 1)))
Usig=matrix(numpy.matlib.zeros((16, 1)))
UF  =matrix(numpy.matlib.zeros((16, 1)))
UFvp=matrix(numpy.matlib.zeros((20, 1)))
Us  =matrix(numpy.matlib.zeros((4, 1)))

#Gauss variables
#d vectors are [xx yx xy yy zz]
F       = matrix(numpy.matlib.zeros((4,1)))#Tensor form
Ft      = matrix(numpy.matlib.zeros((3,3)))#Tensor form
Fd      = matrix(numpy.matlib.zeros((4,1)))#Vector form, FZZ IS 1!!
Ft_inv  = matrix(numpy.matlib.zeros((3,3)))#Tensor form
Fd_inv  = matrix(numpy.matlib.zeros((5,1)))#Tensor form
Fet_inv = matrix(numpy.matlib.zeros((3,3)))#Tensor form
Fet     = matrix(numpy.matlib.zeros((3,3)))#Tensor form

#Formulation 2
Fvp  = matrix(numpy.matlib.zeros((4,1)))#Tensor form
Fvpt = matrix(numpy.matlib.zeros((3,3)))#Tensor form
Fvpd = matrix(numpy.matlib.zeros((5,1)))#Tensor form

S=matrix(numpy.matlib.zeros((4, 1)))    #Nodal internal variable
v=matrix(numpy.matlib.zeros((2, 1)))    #Gauss Point velocity

#Shape Functions
Ns=matrix(numpy.matlib.zeros((1, 4)))
Nv=matrix(numpy.matlib.zeros((2, 8)))
NsigF=matrix(numpy.matlib.zeros((4, 16)))
NFvp=matrix(numpy.matlib.zeros((5, 20)))

Dvp=matrix(numpy.matlib.zeros((4, 1)))      #From plastic deformation direction
DM =matrix(numpy.matlib.zeros((5, 5)))   #4.24
temp2x2=matrix(numpy.matlib.zeros((2, 2)))

#Derivatives
dHxy=matrix(numpy.matlib.zeros((2, 4)))
Bs=matrix(numpy.matlib.zeros((2, 4)))
Bv=matrix(numpy.matlib.zeros((4, 8)))
BsigF=numpy.zeros((4,16,2))
BFvp=numpy.zeros((5,20,2))

temp4x16=matrix(numpy.matlib.zeros((4, 16)))
temp4x1=matrix(numpy.matlib.zeros((4, 1)))
temp5x1=matrix(numpy.matlib.zeros((5, 1)))
temp4x8=matrix(numpy.matlib.zeros((4, 8))) #BL*NFUF
temp5x16=matrix(numpy.matlib.zeros((5, 20)))

temp4x2=matrix(numpy.matlib.zeros((4, 2)))
temp20x2=matrix(numpy.matlib.zeros((20, 2))) #For 4.39
temp16x2=matrix(numpy.matlib.zeros((16, 2))) #For 4.39

#(5,20,2) x (20)
BUFvp=matrix(numpy.matlib.zeros((5, 2)))

B4i=numpy.zeros((4,4,2))
B5i=numpy.zeros((5,5,2))
#(4,16,2)
#print(BsigF[0])
LM =matrix(numpy.matlib.zeros((4, 4)))
LM5=matrix(numpy.matlib.zeros((5, 5)))
#Only returns derivatives respectives to F and Fvp
dEdU=[matrix(numpy.matlib.zeros((4, 16))),matrix(numpy.matlib.zeros((4, 20)))]

K=matrix(numpy.matlib.zeros((44, 44)))

R    =[  matrix(numpy.matlib.zeros(( 8, 1))),
         matrix(numpy.matlib.zeros((16, 1))),
         matrix(numpy.matlib.zeros((20, 1))),
         matrix(numpy.matlib.zeros(( 4, 1)))]
Rzero=[  matrix(numpy.matlib.zeros(( 8, 1))),
         matrix(numpy.matlib.zeros((16, 1))),
         matrix(numpy.matlib.zeros((20, 1))),
         matrix(numpy.matlib.zeros(( 4, 1)))]
RF  =matrix(numpy.matlib.zeros((16, 1)))
Rsig=matrix(numpy.matlib.zeros((16, 1)))
RFvp=matrix(numpy.matlib.zeros((20, 1)))
Rs  =matrix(numpy.matlib.zeros((4, 1)))
Rv  =matrix(numpy.matlib.zeros((8, 1)))

#Formulation 1
#-----------------------------------
#Symmetric tensors
Ee= matrix(numpy.matlib.zeros((4, 1)))
Eet=matrix(numpy.matlib.zeros((3, 3)))
D= matrix(numpy.matlib.zeros((4, 1))) #Strain Rate 4.4
Dt=matrix(numpy.matlib.zeros((2, 2))) #Strain Rate 4.4


#These are the same but reorganized
dVxy=zeros(4)
L   =matrix(numpy.matlib.zeros((2, 2)))
BL      = numpy.zeros((4,4,8))                      #Eqns 2.33, B.17
temp8x1 = matrix(numpy.matlib.zeros((8, 1))) 


#Stress
sig  =matrix(numpy.matlib.zeros((4, 1)))    #Stress Gauss Points xx, yy,zz, syx, is sbar, WHICH IS SYMMETRIC
sig_d=matrix(numpy.matlib.zeros((4, 1)))    #Deviatoric, is also symmetric
signod=matrix(numpy.matlib.zeros((numnodes, numnodes*4)))
N_d  =matrix(numpy.matlib.zeros((4, 1)))    #Direction of plastic strain

P=matrix(numpy.matlib.zeros((4, 1))) 

#Global matrices
#Uglob=matrix(numpy.matlib.zeros((ndof*numnodes, ndof*numnodes)))
Kglob=matrix(numpy.matlib.zeros((dof, dof))) 
dUglob=zeros(dof)
Uglob=zeros(dof)
Uglobant=zeros(dof)
Rglob=zeros(dof)

class bMatrix: #Block matrix
    
    def __init__(self,i,j):
        m=matrix(numpy.matlib.zeros((i, j)))

    # def show_all(self):
        # print(self.left, self.right)
#Matrix Tangent 
#U derivatives
#
# Kt_U=[[matrix(numpy.matlib.zeros(( var_edof.astype(int)[0], var_edof.astype(int)[0]))), matrix(numpy.matlib.zeros((var_edof.astype(int)[0], var_edof.astype(int)[1])))],
      # [matrix(numpy.matlib.zeros(( var_edof.astype(int)[0], var_edof.astype(int)[2]))), matrix(numpy.matlib.zeros((var_edof.astype(int)[0], var_edof.astype(int)[3])))]]
   
#Matrix Double Entry
#Example in formulation 1:
# ( 8x 8)  (8x16)  (8x20) (8x4)
# (16x 8) (16 x 16) ..  ..
# (20x 8) (20x16) (20x20) (20x4)
# ( 4x 8)

dDdUv   =numpy.zeros((5,5, 8)) #TO MODIFY
dDdUF   =numpy.zeros((5,5,16)) #TO MODIFY
dDdUFvp =numpy.zeros((5,5,20)) #TO MODIFY
dDdUs   =numpy.zeros((5,5, 4)) #TO MODIFY

#From multiplying dDdU x NFvp*
#4.39 to 4.42
temp_dDFvp=[numpy.matlib.zeros((5, 8)),
            numpy.matlib.zeros((5,16)),
            numpy.matlib.zeros((5,20)),
            numpy.matlib.zeros((5, 4))]

Kt=[
     [matrix(numpy.matlib.zeros(( var_edof.astype(int)[0], var_edof.astype(int)[0]))), matrix(numpy.matlib.zeros((var_edof.astype(int)[0], var_edof.astype(int)[1]))),
      matrix(numpy.matlib.zeros(( var_edof.astype(int)[0], var_edof.astype(int)[2]))), matrix(numpy.matlib.zeros((var_edof.astype(int)[0], var_edof.astype(int)[3])))]
     ,
     [matrix(numpy.matlib.zeros(( var_edof.astype(int)[1], var_edof.astype(int)[0]))), matrix(numpy.matlib.zeros((var_edof.astype(int)[1], var_edof.astype(int)[1]))),
      matrix(numpy.matlib.zeros(( var_edof.astype(int)[1], var_edof.astype(int)[2]))), matrix(numpy.matlib.zeros((var_edof.astype(int)[1], var_edof.astype(int)[3])))
     ],
     [matrix(numpy.matlib.zeros(( var_edof.astype(int)[2], var_edof.astype(int)[0]))), matrix(numpy.matlib.zeros((var_edof.astype(int)[2], var_edof.astype(int)[1]))),
      matrix(numpy.matlib.zeros(( var_edof.astype(int)[2], var_edof.astype(int)[2]))), matrix(numpy.matlib.zeros((var_edof.astype(int)[2], var_edof.astype(int)[3])))]
     ,
     [matrix(numpy.matlib.zeros(( var_edof.astype(int)[3], var_edof.astype(int)[0]))), matrix(numpy.matlib.zeros((var_edof.astype(int)[3], var_edof.astype(int)[1]))),
      matrix(numpy.matlib.zeros(( var_edof.astype(int)[3], var_edof.astype(int)[2]))), matrix(numpy.matlib.zeros((var_edof.astype(int)[3], var_edof.astype(int)[3])))
     ]    
    ] 
    
Kzero=[
     [matrix(numpy.matlib.zeros(( var_edof.astype(int)[0], var_edof.astype(int)[0]))), matrix(numpy.matlib.zeros((var_edof.astype(int)[0], var_edof.astype(int)[1]))),
      matrix(numpy.matlib.zeros(( var_edof.astype(int)[0], var_edof.astype(int)[2]))), matrix(numpy.matlib.zeros((var_edof.astype(int)[0], var_edof.astype(int)[3])))]
     ,
     [matrix(numpy.matlib.zeros(( var_edof.astype(int)[1], var_edof.astype(int)[0]))), matrix(numpy.matlib.zeros((var_edof.astype(int)[1], var_edof.astype(int)[1]))),
      matrix(numpy.matlib.zeros(( var_edof.astype(int)[1], var_edof.astype(int)[2]))), matrix(numpy.matlib.zeros((var_edof.astype(int)[1], var_edof.astype(int)[3])))
     ],
     [matrix(numpy.matlib.zeros(( var_edof.astype(int)[2], var_edof.astype(int)[0]))), matrix(numpy.matlib.zeros((var_edof.astype(int)[2], var_edof.astype(int)[1]))),
      matrix(numpy.matlib.zeros(( var_edof.astype(int)[2], var_edof.astype(int)[2]))), matrix(numpy.matlib.zeros((var_edof.astype(int)[2], var_edof.astype(int)[3])))]
     ,
     [matrix(numpy.matlib.zeros(( var_edof.astype(int)[3], var_edof.astype(int)[0]))), matrix(numpy.matlib.zeros((var_edof.astype(int)[3], var_edof.astype(int)[1]))),
      matrix(numpy.matlib.zeros(( var_edof.astype(int)[3], var_edof.astype(int)[2]))), matrix(numpy.matlib.zeros((var_edof.astype(int)[3], var_edof.astype(int)[3])))
     ]    
    ] 
        
dgdU=[  matrix(numpy.matlib.zeros((1, 8))),
        matrix(numpy.matlib.zeros((1,16))),
        matrix(numpy.matlib.zeros((1,20))),
        matrix(numpy.matlib.zeros((1, 4)))]

print("Kt_0",len((Kt[1][0])),len((Kt[1][0]).transpose()))


ck = ey/ ((1. + nu)*(1. - 2. * nu))
c=matrix(numpy.matlib.zeros((4, 4)))
#PLAIN STRAIN
c[0,0]=c[1,1]=c[2,2]=ck*(1.-nu)                 
c[0,1]=c[1,0]=c[0,2]=c[2,0]=c[2,1]=c[1,2]=ck*nu
c[3,3]=ck*(.5 - nu)
print ("C matrix")
print(c)
#Set G and H matrices
#2.32
G=matrix([[1,0,0,0],[0,0,0,1],[0,0,0,1],[0,1,0,0]]) 
H=matrix([[1,0,0,0],[0,1,0,0],[0,0,0,0],[0,0,0,0]]) 


#---------------------------------------------------
#Before start, assign Uglob values
#---------------------------------------------------
#To MODIFY put in a function     

# for n in range(numnodes):
    # #Velocities 
    # iF=ndof*n
    # Uglob[iF  ]=    Uglob[iF +1] = 0. #Velocities 
    # #Initial deformation gradients as identity??
    # Uglob[iF+2]=Uglob[iF+5]=1
    # Uglob[iF+3]=Uglob[iF+4]=0#xy and yx        	
# print ("Initial Uglob", Uglob)

#TODO: SHORTEN THIS
for inode in range(numnodes):
    #print("node",inode)   
    #Deformation gradient F
    vrowinc=0
    varinc=0
    for nvar in range(numvars):
        print("nvar, vrowinc",nvar, vrowinc)
        idof=int(vrowinc+var_dim[nvar]*inode)
        for i in range ( var_dim [ nvar ] ):
            Uglob[idof] = Uinit[varinc+i] 
            idof+=1
        varinc+=var_dim[nvar]
        vrowinc+=numnodes*var_dim[nvar]      

for n in range(size(boundarynode)):
    inode=boundarynode[n]
    print("node",inode)   
    #Deformation gradient F
    vrowinc=0
    varinc=0
    for nvar in range(numvars):
        print("nvar, vrowinc",nvar, vrowinc)
        idof=int(vrowinc+var_dim[nvar]*inode)
        for i in range ( var_dim [ nvar ] ):
            if is_bcnode_byvar[n,nvar]:
                print ("BC idof, varinc+i, val",idof,varinc+i,node_bc[ n, varinc+i ])
                if solver == 2:
                    Uglob[idof] = node_bc[ n, varinc+i ]           #F INCREMENT (dF) IS NULL!!!!!
            else:
                print ("INITIAL idof, varinc+i,val",idof,varinc+i,Uinit[varinc+i] )
                Uglob[idof] = Uinit[varinc+i] 
            idof+=1
        varinc+=var_dim[nvar]
        vrowinc+=numnodes*var_dim[nvar]                           
                
print ("Initial Uglob (with bcs)",Uglob)
#-------------------------------------------
it=0

def calculate_shape_matrices (rg, sg,X2):
    #rg=gauss[ig]
    #sg=gauss[jg]

    #Numerated as in Bathe
    Ns  =0.25*matrix([(1+sg)*(1+rg),(1-rg)*(1+sg),(1-sg)*(1-rg),(1-sg)*(1+rg)])   
    dHrs=matrix([[(1+sg),-(1+sg),-(1-sg),(1-sg)], [(1+rg),(1-rg),-(1-rg),-(1+rg)] ])
    #Numerated as in deal.ii
    #dHrs=matrix([[-(1-s),(1-s),-(1+s),(1+s)], [-(1-r),-(1+r),(1-r),(1+r)] ])        
    dHrs/=4.
    J=dHrs*X2
    dHxy=linalg.inv(J)*dHrs
    detJ=linalg.det(J)
    #Calculate shape functions
    #Bs=J-1 dHrs(B.13)
    Bs=dHxy
    for k in range(4):
        #shape functions
        Nv[0,2*k  ]=Nv[1,2*k+1]=Ns[0,k]
        for j in range(4):
            NsigF[j,4*k+j]=Ns[0,k]

        #derivatives Bv (B.14)
        Bv[0,2*k  ]=dHxy[0,k]
        Bv[1,2*k  ]=dHxy[1,k]
        Bv[2,2*k+1]=dHxy[0,k]
        Bv[3,2*k+1]=dHxy[1,k]
    
    ################            
    ## 
    for i in range(4):
        for l in range(4):
            for m in range(4):  
                for n in range(2):
                    if (l==m):
                        B4i[l,m,n]=Bs[n,i]
                    else:
                        B4i[l,m,n]=0.              
        for l in range(4):
            for m in range(4):  
                for n in range(2):
                    BsigF[l,4*i+m,n]=B4i[l,m,n]
    
    ###################
    ## Ec. D.20 p177 
    if form==2:
        for i in range(4):
            for l in range(5):
                for m in range(5):  
                    for n in range(2):
                        if (l==m):
                            B5i[l,m,n]=Bs[n,i]
                        else:
                            B5i[l,m,n]=0.              
            for l in range(5):
                for m in range(5):  
                    for n in range(2):
                        BFvp[l,4*i+m,n]=B5i[l,m,n]   

    return Nv,NsigF,Ns                    

## ------------------------------------------
## Newton Rhapson Loop
## ------------------------------------------
if plastic == 0:
    UFvp=Fvpd=[1.,0.,0.,1.,1.]
print ("**************************************** BEGIN LOOP *************************************************")
while (it < numit):
    print ("Iteration: ",it, "from ", numit)
    #Clean Global Matrices for assembly
    #print ("Kglob",Kglob)
    for idof in range(dof):
        Rglob [idof] = 0.
        for jdof in range(dof):
            Kglob[idof,jdof] = 0.  
#ELEMENT LOOP  ----------------
    for e in range (numel): # TO MODIFY 
        #Obtain Ve from global
        Kt=Kzero
        R=Rzero
        Kel=0.
        for n in range(4):
            X2[n]=node[elnodes.astype(int)[e][n]]
        #print ("Element ", e)
        #print ("Element Nodes")
        #print (X2)
        
        dHrs=matrix([[(1+0.),-(1+0.),-(1-0.),(1-0.)], [(1+0.),(1-0.),-(1-0.),-(1+0.)] ])
        J=dHrs*X2
        lx=J[0,0]/2.
        ly=J[1,1]/2.
        #print("J,lx,ly",J,lx,ly)
        for ig in range(2):
            for jg in range(2):
                rg=gauss[ig]
                sg=gauss[jg]
                detJ=0.
                # Nv,NsigF,Ns=calculate_shape_matrices (rg, sg,X2)
                # print ("dHxy",dHxy)
                #Numerated as in Bathe
                Ns  =0.25*matrix([(1+sg)*(1+rg),(1-rg)*(1+sg),(1-sg)*(1-rg),(1-sg)*(1+rg)])   
                dHrs=matrix([[(1+sg),-(1+sg),-(1-sg),(1-sg)], [(1+rg),(1-rg),-(1-rg),-(1+rg)] ])
                #Numerated as in deal.ii
                #dHrs=matrix([[-(1-s),(1-s),-(1+s),(1+s)], [-(1-r),-(1+r),(1-r),(1+r)] ])        
                dHrs/=4.
                J=dHrs*X2
                dHxy=linalg.inv(J)*dHrs
                detJ=linalg.det(J)
                #Calculate shape functions
                #Bs=J-1 dHrs(B.13)
                Bs=dHxy
                for k in range(4):
                    #shape functions
                    Nv[0,2*k  ]=Nv[1,2*k+1]=Ns[0,k]
                    for j in range(4):
                        NsigF[j,4*k+j]=Ns[0,k]

                    #derivatives Bv (B.14)
                    Bv[0,2*k  ]=dHxy[0,k]
                    Bv[1,2*k  ]=dHxy[1,k]
                    Bv[2,2*k+1]=dHxy[0,k]
                    Bv[3,2*k+1]=dHxy[1,k]
                
                ################            
                ## 
                for i in range(4):
                    for l in range(4):
                        for m in range(4):  
                            for n in range(2):
                                if (l==m):
                                    B4i[l,m,n]=Bs[n,i]
                                else:
                                    B4i[l,m,n]=0.              
                    for l in range(4):
                        for m in range(4):  
                            for n in range(2):
                                BsigF[l,4*i+m,n]=B4i[l,m,n]
                
                ###################
                ## Ec. D.20 p177 
                if form==2:
                    for i in range(4):
                        for l in range(5):
                            for m in range(5):  
                                for n in range(2):
                                    if (l==m):
                                        B5i[l,m,n]=Bs[n,i]
                                    else:
                                        B5i[l,m,n]=0.              
                        for l in range(5):
                            for m in range(5):  
                                for n in range(2):
                                    BFvp[l,4*i+m,n]=B5i[l,m,n]            
                                
                #Interpolate velocity
                #INCREMENT GLOBAL VELOCITY FROM INCREMENTS!!!
                #CHANGE F TO ASSEMBLE IN THE SAME PLACE FOR BOTH FORMS
                # THIS IS THE CRITERIA IN WHICH VARS ARE INBLOCK PER NODE
                # juf=0
                # uvf=0
                # for n in range (4):
                    # d=elnodes.astype(int)[e][n]
                    # for i in range (var_dim[0]):    #Velocity is var 0
                        # print("UV loc glob ",i,ndof*d+i)
                        # UV[i,0]=Uglob[ndof*d+i]
                    # uvf+=var_dim[0]
                    # for j in range (var_dim[1]):
                        # #print("J",j)
                        # if (form==1):
                            # Usig[j+juf,0]=Uglob[ndof*d+var_dim[0]+j]
                            # UF  [j+juf,0]=Uglob[ndof*d+6+j]
                        # else: #Fig 4.1, Z is not translated to Fvpt
                            # UF  [j+juf,0]=Uglob[ndof*d+var_dim[0]+j]
                            # #print("UF(j,coord)",j,ndof*d+6+j)
                    # juf+=var_dim[1]
                    
                    # if plastic:
                        # for j in range (var_dim[2]):
                            # UFvp[j,0]=Uglob[ndof*d+var_dim[0]+var_dim[1]+j]
               
                #IN CASE GLOBAL NODES ARE STORED IN BLOCK
                vrowinc=0
                #Assembly Matrix
                for vrow in range(numvars): #Variables
                    #print("vrow",vrow)
                    ir=0
                    imax=int(var_dim[vrow])
                    for n in range (4): #Nodes
                        for i in range(imax): 
                            d=elnodes.astype(int)[e][n]
                            #print("ir glob",ir, vrowinc+var_dim[vrow]*d+i)
                            vnrow[ir]=vrowinc+var_dim[vrow]*d+i
                            ir=ir+1
                                
                        # print("vnrow",vnrow.astype(int)) 
                    
                    if   vrow == 0:
                        for row in range(4*imax):
                            UV[row,0]=Uglob[int(vnrow[row])]
                    elif vrow == 1:
                        for row in range(4*imax):
                            UF[row,0]=Uglob[int(vnrow[row])]                                
                    vrowinc+=numnodes*var_dim[vrow]

                            
                #print("UF",UF)

                v  =Nv*UV #[2x8 x (8x1)]
                s  =float(Ns*Us)
                F  =NsigF*UF #[(4x16)*(16x1) =(4x1)]
                #print("F",F)
                #Ft=
                if (form==1):
                    sig=NsigF*Usig
                else:
                    if plastic:
                        Fvp=NFvp*UFvp
                
                 
                #Galerkin strain integration
                #Calculate deformation gradient Fij, for that
                #Calculate Velocity gradient Lij
                #Lij=dvi/dxj(2.4) 
                #According to B.11 d(phi)/dxj=J-1(ij) dphi/drj = Bvjk Vk
                
                #Formulation 1
                #Calculate La (ij) 2.9 
                #Laij=F(-1)ki F(-1)kl Le lj
                
                #Calculate Leij = Lij - D(th)ij - D(vp) ij (2.10-2.12)
                #Calculate Almansi deformation gradient E (A.5)
                #Calculate 
                #Ea ij= 1/2(Lki F(-1)lk F(-1)LJ +F(-1)ki F(-1)KL Llj )
                dVxy=Bv*UV #(4x8)*(8x1)=(4x1) (vx,x vx,y vy,x vy,y)T 
                L[0,0]=dVxy[0]
                L[0,1]=dVxy[1]
                L[1,0]=dVxy[2]
                L[1,0]=dVxy[3]
                
                Dt=0.5*(L+L.transpose())
                D[0]=Dt[0,0];D[0]=Dt[0,0];
                
                #Stabilization factor tau 2.26
                #tau=beta*he/(2|v|)
                #See beta estability paramter
                #LM (2.33 p23)
                #Attention: SE DIFFERENCES WITH L_ in 2.28
                LM[0,0]=LM[2,2]=dVxy[0]
                LM[0,1]=LM[2,3]=dVxy[1]
                LM[1,0]=LM[3,2]=dVxy[2]
                LM[1,1]=LM[3,3]=dVxy[3]
                
                for i in range(4):
                    for j in range(4):               
                        LM5[i,j]=LM[i,j]
                        
                #BL interpolators BLijk (4,4,8) (B.17 p165)
                for k in range(8):
                    BL[0,0,k]=BL[2,2,k]=Bv[0,k]
                    BL[0,1,k]=BL[2,3,k]=Bv[2,k]
                    BL[1,0,k]=BL[3,2,k]=Bv[1,k]
                    BL[1,1,k]=BL[3,3,k]=Bv[3,k]
                    BL[0,2,k]=BL[0,3,k]=BL[1,2,k]=BL[1,3,k]=BL[2,0,k]=BL[2,1,k]=BL[3,0,k]=BL[3,1,k]=0.
                #Set tractions tP (2.34)
                
                
                    
                w=1. #TO MODIFY
                
                beta=1.
                #print("u2+v2",v[0]*v[0]+v[1]*[1])
                vnorm=sqrt(v[0]*v[0]+v[1]*v[1])
                #he=(lx+ly)/2.
                he=(lx*v[0]+ly*v[1])/vnorm #HUGHES (APPROX)
                tau=float(beta*he/(2.*vnorm))

                Ft[0,0]=F[0]
                Ft[0,1]=F[1]
                Ft[1,0]=F[2]
                Ft[1,1]=F[3]
                Ft[2,2]=1.      #Fzz, plain strain 
                
                #F is [xx xy yx yy zz] , (4.21) in form 2
                Fd[0]=F[0]
                Fd[1]=F[2] #yx
                Fd[2]=F[1] #xy
                Fd[3]=F[3]
                                
                #Calculate stabilization parameter
                tau=1.
                ################ STRESSES**********
                #FORMULATION TRUE EQUILIBRIUM, EQ 4.2
                #####################################
                
                #FORM 1
                #Calculate sigma
                #2.31 Comes first from 2.2 (Then 2.31)
                #Pij=vk d(sig ij)/xk - dvi/dxk sig(kj) + dvk/dxk sig ij
                #Calculate Piola Kirchoff Pi (2.31) Gij Cjk˙¯-Gij LM (jk) sig(k) +
                #Attention double contraction
                #P=G*c*E-G*LM*sig+(LM[0,0]+LM[1,1])*G*sig
                
                #FORM 2
                if plastic == 0:
                    Fet=Ft
                
                Fet_inv=linalg.inv(Fet) 
                Eet=0.5*(identity(3)-Fet_inv.transpose()*Fet_inv)
                #4 x 1 vector arrangement, for constitutive equation, Eq. 4.20
                Ee[0]=Eet[0,0];Ee[1]=Eet[1,1];Ee[2]=Eet[2,2];
                Ee[3]=Eet[0,1]; #OR 1,0, being Eet symmetric
                visc=1.
                sig=c*Ee+2*visc*D
                #print ("sig",sig)
                
                #E if were total Is the same as Almansi?? (NONLINEAR CONTINUA EQ. 4.3)
                
                #AND THEN 4.20
                #From 2.27 Plane Strain Symmetric tensors are defined as 
                #t=[txx tyy tzz tyz]
                pi=1./3.*(sig[0,0]+sig[1,0]+sig[2,0])
                #print("pi",pi)
                for i in range(3): #Only daigonal is modified
                    sig_d[i,0]=sig[i,0]-pi #comps are [x y z yz]
                #print ("sigd",sig_d[i][0])   
                for k in range(4):
                    sig_eq=sqrt(1.5*(sig_d[k,0]))
                #print ("sig_eq",sig_eq)
                #*** STRAINS
                #Equivalent strain rate
                mat_A=mat_A0*math.exp(-mat_Q/mat_R)
                #Equivalent viscoplastic strain rate
                epsr_vp_eq=0.
                if (s!=0):
                    epsr_vp_eq=mat_A*(sinh(mat_psi*sig_eq/s))**(1./mat_m)
                
                #print (epsr_vp_eq)
                #Evaluate g function 2.58/4.48
                g_sigs=1.
                #g_sigs=h0*|(1-s/s*)|*sign(1-s/s*)A(sinh())
                #With s*
                #s*=s~(edot~_vp/A)^n
                #Direction of plastic Strain: 2.16
                N_d=sqrt(3./2.)*sig_d/sig_eq
                Dvp=sqrt(3./2.)*epsr_vp_eq*N_d
                #print ("Nd",N_d)
                #With evp = f(sigma,s) (Sec. 2.4.2)
                #eps_vp=A sinh (psi (sigma/s))^(1/m) (2.57 p27)
                #Calculate Rate of Def Tensor D (2.13, 2.14)
                #D(th)ij=alpha vk dT/dxk deltaij 2.13
                #IN THIS EXAMPLE THIS DTH 0
                #D(vp)ij=sqrt(3/2) e. vp Nij  2.14
                #Nij Direction of plastic flow
                #Assemble DM matrix 4.24 p92
                if form==2:
                    #for i in range(2):
                    #    for j in range(2):
                    #Dvp comes from N, which comes from 
                    #Dvp is [x y z yz]
                    DM[0,0]=DM[2,2]=Dvp[0]
                    DM[1,1]=DM[3,3]=Dvp[1]
                    DM[1,0]=DM[3,2]=DM[0,1]=DM[2,3]=Dvp[3]
                    DM[3,3]=Dvp[2]
                
                wJ=w*detJ
                
                
                #Elastic part of F 4.6
                #ATENTION F~ is not in the same order of NODAL variable 
                if plastic:
                    if (it==0):
                        #Ft=identity(3)
                        #NOT USE!!! Fd=[[1,0,0,1]]
                        Fvpt=identity(3)
                        #print(Ft)
                        Fvpd[0]=1.
                        Fvpd[1]=0. #yx
                        Fvpd[2]=0. #xy
                        Fvpd[3]=1. 
                        Fvpd[4]=1.
                    else:
                        #Fvp is [xx xy yx yy zz]
                        Fvpd[0]=Fvp[0]
                        Fvpd[1]=Fvp[2] #yx
                        Fvpd[2]=Fvp[1] #xy
                        Fvpd[3]=Fvp[3]  
                        Fvpd[4]=Fvp[4]
                    
                
                #print ("Fvpd",Fvpd)
                #print("Fd",Fd[0,0])
                #Arguments passed are ~ vectors
                dEdU=deriv.calc_dEdU(Fd,Fvpd,NsigF,NFvp)
                
                # *****************************************************************
                #RESIDUALS ******************* 2.26 to 2.39 *****************************
                R[0]  =Bv.transpose()*P #Remains the summ of particular gauss points
                #Rsig[16x1] (4 per node)
                #Construct vk Bsig mik
                for m in range(4):
                    for i in range(16):
                            temp4x16[m,i]=0.
                for m in range(4):
                    for i in range(16):
                        for k in range(2):
                            temp4x16[m,i]=temp4x16[m,i]+BsigF[m,i,k]*v[k,0]
                
                if form==2:
                    for m in range(5):
                        for i in range(20):
                            for k in range(2):
                                temp5x16[m,i]=temp5x16[m,i]+BFvp[m,i,k]*v[k,0]
                            
                if (form==1):
                    Rsig=(NsigF+temp4x16*tau).transpose()*(temp4x16*Usig-c*Ee)*wJ
                else: #4.29
                    if plastic == 1:
                        R[2]=(NFvp+temp5x16*tau).transpose()*(temp5x16*UFvp-DM*NFvp*UFvp)*wJ
                
                R[1]   =(NsigF+temp4x16*tau).transpose()*(temp4x16*UF-LM*NsigF*UF)*wJ
                R[3]    =(Ns+tau*v.transpose()*Bs).transpose()*(v.transpose()*Bs*Us-g_sigs)*wJ
                
                
                #R Assembly            
                #TANGENT MATRIX   
                #PAGES 25 y 94
                #dRdUn=Bv.transpose()*( G*c*dEdU -G*)
                # for j,n in range(8,8):
                    # for i in range(8,8):
                    # Kt[0][0]=   Kt[0][0]+B[i][j]*(
                                # #*G[i][p]*C[p][k]
                                # -G[i][l]*BL[l][k][n]*sigma[k]
                                # +G[i][k]*sigma[k]*H[m][l]*BL[m][l][n]*
                                # wJ)   
                #temp8x1=0
                #*** ONLY FOR FORMULATION 2
                # for l,m,n in range(4,4,8):
                    # temp8x1[l]=temp8x1[l]+BL[m][l][n]
                
                # #print("temp",temp8x1[0])  
                # Kt[0][0]=   Kt[0][0]+Bv.transpose()*(
                            # #G[i][p]*C[p][k]
                            # #-G[i][l]*BL[l][k][n]*sigma[k]
                            # +G*sig*temp8x1.transpose()*
                            # wJ
                #dRv/dUv
                Kt[0][0]=   Kt[0][0]+Bv.transpose()*visc*Bv*wJ
                
                #dRv/dUF
                #(8x16)
                Kt[0][1]=   Kt[0][1]+Bv.transpose()*c*dEdU[0]*wJ
                #dRv/dUF (8x20)
                if plastic:
                    Kt[0][2]=   Kt[0][2]+Bv.transpose()*c*dEdU[1]*wJ
                #Kt(0,3) dRv/dUs=0
                
                #----------------------------------------------------
                #for i in range(4):
                #Kt.clear()
                #F derivatives (2.49 to 2.52)
                if (form==1):
                    find=int(2)
                else:
                    find=int(1)
                    
                for m in range(4):
                    for j in range(16):
                        for k in range(2):
                            temp4x2[m,k]=temp4x2[m,k]+BsigF[m,j,k]*UF[j,0]

                for m in range(4):
                    for l in range (4):
                        for p in range(8):
                            temp4x8[m,p]=BL[m,l,p]+F[l,0]            
                #dRF/dUv 4.35
                #(16x8)
                #print("Nv",Nv)
                
                for m in range(4):
                    for i in range(16):
                        for n in range(2):
                            temp1=temp2=0. #increase with each m
                            
                            for j in range(16):
                                for k in range(2):
                                    temp1=temp1+2*BsigF[m,j,k]*v[k]*UF[j,0]
                            
                            temp16x2[i,n]=  temp20x2[i,n]+(
                                            BFvp[m,i,n]*
                                            (temp1-temp4x1[m,0])
                                            )
                                            
                Kt[find][0]   =Kt[find][0]+( 
                                 (NsigF.transpose()*temp4x2*Nv) #16x4*4x2*2x8
                                +tau*temp16x2*Nv #temp_in  * N_np
                                -(NsigF+float(tau)*temp4x16).transpose()*temp4x8
                                )*wJ
                #dRFdUF  4.36
                Kt[find][find]= Kt[find][find]+(
                                (NsigF+temp4x16*tau).transpose()*
                                (temp4x16-LM*NsigF)
                                )*wJ
                                
                #dRFdUF=dRFdUF=0.  4.37 & 4.38
                
                if plastic:
                    #---------------------------------------------------
                    ##Viscoplastic derivatives
                    #Kt(2,0)=dFvp/dUV 4.39
                    temp5x1=LM5*NFvp*UFvp
                    temp1=0
                    #print ("temp4x1,temp1",temp5x1,temp1)
                    
                    for m in range(5):
                        for i in range(20):
                            for n in range(2):
                                temp1=temp2=0. #increase with each m
                                
                                for j in range(20):
                                    for k in range(2):
                                        temp1=temp1+2*BFvp[m,j,k]*v[k]*UFvp[j,0]
                                
                                temp20x2[i,n]=  temp20x2[i,n]+(
                                                BFvp[m,i,n]*
                                                (temp1-temp5x1[m,0])
                                                )
                    
                    for m in range(5):
                        for j in range (20):
                            for k in range (2):
                                BUFvp[m,k]=BFvp[m,j,k]*UFvp[j,0]
                    #print("size",len((NFvp+float(tau)*temp5x16)))
                    #dRFvp/DUv
                
                    Kt[2][0]=Kt[2][0]+(
                                NFvp.transpose()*BUFvp*Nv+
                                 tau*temp20x2*Nv #temp_in  * N_np
                                -(NFvp+float(tau)*temp5x16).transpose()*temp_dDFvp[0]
                                )*wJ
                    
                    #dRFvp/dUF 4.40
                    Kt[2][1]=Kt[2][1]-(
                                 (NFvp+float(tau)*temp5x16).transpose()*temp_dDFvp[1]
                                )*wJ
                    
                    #dRFvp/dUFvp 4.41
                    Kt[2][2]=Kt[2][2]+(
                                 (NFvp+float(tau)*temp5x16).transpose()*(
                                 temp5x16-DM*NFvp-temp_dDFvp[2])
                                )*wJ

                    #dRFvp/dUs 4.42
                    Kt[2][3]=Kt[2][3]-(
                                 (NFvp+float(tau)*temp5x16).transpose()*(
                                 temp_dDFvp[3])
                                )*wJ
                                
                    #---------------------------------------------------
                    #S derivatives -- COMMON TO BOTH FORMULATIONS
                    #2.53 and 4.xx ATENTION IN CHAPTER 4 k and p indices are wrong, see 2.53
                    #dRs/dUv
                    # print("test",tau*Bs.transpose()*Nv*(2.*v*Bs*Us * g_sigs) 
                    #print("test",float(2.0*(Bs*Us).transpose()*v))
                    Kt[3][0]=Kt[3][0]+(
                              Ns.transpose()*(Bs*Us).transpose()*Nv +
                              tau*Bs.transpose()*Nv*float(2.*(Bs*Us).transpose()*v * g_sigs)
                              )*wJ 
                    
                    #4.44 dRs/dF 4.44 (4x16)
                    Kt[3][1]=Kt[3][1]+(
                              (Ns+tau*v.transpose()*Bs).transpose()*
                              (-dgdU[1])*
                              wJ) 
                    
                    #dRs/dFvp 4.45
                    #4x20
                    Kt[3][2]=Kt[3][2]+(
                            (Ns.transpose()+tau*Bs.transpose()*v)*(-1.)*dgdU[2]
                            )*wJ
                    #dRs/dS 4.46    
                    Kt[3][3]=Kt[3][3]+(
                              (Ns+tau*v.transpose()*Bs).transpose()*
                              (v.transpose()*Bs-dgdU[3])*
                              wJ) 
                #END IF PLASTIC ******************************************
                
        #print ("Nv",Nv)
        # #Element dimension and DOF PER VARIABLE! 
        # var_dim =[2,4,4,1]
        # var_edof=zeros(4)
        # if form==2:
            # var_dim =[2,4,5,1]
        # for i in range(4):
            # var_edof[i]=4*var_dim[i]  
            
        #Assuming alternated variables
        #From https://www.dealii.org/current/doxygen/deal.II/step_20.html
        #However, then things become different. As mentioned in the introduction, we want to subdivide the matrix into blocks corresponding 
        #to the two different kinds of variables, velocity and pressure. To this end, we first have to make sure that the indices corresponding 
        #to velocities and pressures are not intermingled: 
        #First all velocity degrees of freedom, then all pressure DoFs. This way, the global matrix separates nicely into a 2×2 system. 
        #To achieve this, we have to renumber degrees of freedom base on their vector component, an operation that conveniently is already implemented:
        #Example 9 nodes 4 elements (bad node numbering)
        # Var block numbering
        # [12 13] [14 15] [16 17]
        # [6 7, 30..33]  [8 9,34..37] [10 11, 38.41]
        # [0 1, 18..21]  [2 3,22..25] [4 5, 26..29] 
        
        #Var 0
        # Nodes 4 3 0 1
        # vn=[8 9 6 7 0 1 2 3]
        # vn=[34..37 30.33 18..21]
            
        vrowinc=0
        #Assembly Matrix
        for vrow in range(numvars): #Variables
            #print("vrow",vrow)
            ir=0
            imax=int(var_dim[vrow])
            for n in range (4): #Nodes
                for i in range(imax): 
                    d=elnodes.astype(int)[e][n]
                    #print("vrowinc,d,a+b",vrowinc,d,vrowinc+var_dim[vrow]*d+i)
                    vnrow[ir]=vrowinc+var_dim[vrow]*d+i
                    ir=ir+1
            
            vcolinc=0        
            for vcol in range(numvars): #Variables
                #print("vcol",vcol)
                jmax=int(var_dim[vcol])
                #print("imax, jmax",imax,jmax)
                #Store vn vectors
                ic=0
                for n in range (4): #Nodes 
                    for j in range(jmax):
                        d=elnodes.astype(int)[e][n]
                        #print("vcolinc",vcolinc)
                        vncol[ic]=vcolinc+var_dim[vcol]*d+j
                        ic=ic+1
                            
                        
                #print("vnrow",vnrow.astype(int))            
                #print("vncol",vncol.astype(int))
                for row in range(4*imax):
                    Rglob[vnrow.astype(int)[row]]+=R[vrow][row]
                    for col in range(4*jmax):
                        Kglob[vnrow.astype(int)[row],vncol.astype(int)[col]]=  Kglob[vnrow.astype(int)[row],vncol.astype(int)[col]]+(
                                                                              Kt[vrow][vcol][row,col])
                vcolinc+=numnodes*var_dim[vcol]
            
            
            vrowinc+=numnodes*var_dim[vrow]
        #print (K)

    ##---------------------------------------------------------------
    ##Boundary conditions
    ##---------------------------------------------------------------
    #Velocity DOFs 
    #In this example velocities are known
    #AT INLET(left nodes):
    # F=I , sigma = 0
    if solver == 1: #NONZERO VALUES
        for n in range(size(boundarynode)):
            inode=boundarynode[n]
            idof = ndof * inode 
            for nvar in range(numvars):
                for i in range ( var_dim [ nvar ] ):
                    #print("idof",idof)
                    if is_bcnode_byvar[n,nvar]:
                        for j in range(dof):
                            Rglob[ j ] = Rglob[ j ] - Kglob[j,int(idof)] * node_bc[ n, i ] #dU=0, U=1(idof)
                    idof+=1 #Increase var always although var is not bc
        
        for n in range(size(boundarynode)):
            inode=boundarynode[n]
            for i in range ( var_dim [ 0 ] ):
                idof = var_dim[0] * inode + i   
                #print("i, idof",i, idof)
                Rglob[ int(idof) ] = node_bc[ n, i ] #dU=0, U=1(idof)    


    # ORIGINAL BC APPLICATION; IN CASE OF BLOCK GLOBAL DOFS ER EACH NODE
    # for n in range(size(boundarynode)):
        # inode=boundarynode[n]
        # idof = int(ndof * inode)
        # print("node",inode)   
        # #Deformation gradient F
        # for nvar in range(numvars):
            # for i in range ( var_dim [ nvar ] ):
                # if is_bcnode_byvar[n,nvar]:
                    # print ("idof",idof)
                    # for j in range(dof):
                        # Kglob[ idof , j ] = 0
                        # #Rglob[ j ] -= Kglob[j,idof] * 0 #dU=0 in NEWTON RAPHSON, U=1(idof) IN PICARD
                        # Kglob[ j ,idof ] = 0
                    
                    # Kglob[idof,idof] = 1
                    # if solver == 2:
                        # Rglob[idof  ] = 0           #F INCREMENT (dF) IS NULL!!!!!   
                # idof+=1
    
    #CASE FR SECUENTIAL DOF (FIRST DOFS VELOCITY,THEN F OR SIGMA; etc
    for n in range(size(boundarynode)):
        inode=boundarynode[n]
        #print("node",inode)   
        #Deformation gradient F
        vrowinc=0
        for nvar in range(numvars):
            #print("nvar, vrowinc",nvar, vrowinc)
            idof=int(vrowinc+var_dim[nvar]*inode)
            for i in range ( var_dim [ nvar ] ):
                if is_bcnode_byvar[n,nvar]:
                    #print ("idof",idof)
                    for j in range(dof):
                        Kglob[ idof , j ] = 0
                        #Rglob[ j ] -= Kglob[j,idof] * 0 #dU=0 in NEWTON RAPHSON, U=1(idof) IN PICARD
                        Kglob[ j ,idof ] = 0
                    
                    Kglob[idof,idof] = 1
                    if solver == 2:
                        Rglob[idof  ] = 0           #F INCREMENT (dF) IS NULL!!!!!   
                idof+=1
        
            vrowinc+=numnodes*var_dim[nvar]                           
                
    #print("KGLOB\n")
    # for i in range (dof):
        # for j in range (dof):
            # print(Kglob[i,j], end = " ")
        # print("\n")   
        # print("\n")
    #print("Rglob",Rglob)
    

#print (K)
    if   solver == 1:
        for i in range(dof):
            Uglobant[i]=Uglob[i]
        Uglob=linalg.solve(Kglob, Rglob)
    elif solver == 2:
        dUglob=linalg.solve(Kglob, Rglob)

        for i in range (dof):
            Uglob[i]=Uglob[i]+dUglob[i]
        
    max=0.
    if solver ==1:
        for i in range (dof):
            if abs(Uglobant[i]-Uglob[i])>max:
                max=abs(Uglobant[i]-Uglob[i])
    elif solver ==2:
        for i in range (dof):
            if abs(dUglob[i])>max:
                max=abs(dUglob[i])
    
    print("max dU=",max)
    #print("it %d, dUglob",it, dUglob)  
    #print("Uglob", Uglob)
    
    
    #TOTAL BOUNDARY CONDITIONS FOR UF and s calculations
    # dnode=(nex+1)    
    # for dy in range(ney+1): 
        # inode=dy*dnode
        # #Deformation gradient F
        # idof = var_dim[0] * inode 
        # print ("inode, idof",inode, idof)
        # Uglob[ idof     ] = Uglob[ idof + 3 ] = 1
        # Uglob[ idof + 1 ] = Uglob[ idof + 2 ] = 0

        
        # #Sigma is zero, Internal variable s ,      
        # if numvars == 2:
            # idofs = idof + var_dim[0]
            # Uglob[ idofs ] = mat_s0
 
    
    it+=1
    


    
#print ("Results")
#print("Uglob", Uglob)

varname=["Veloc","DefGrad_F","Fvp","s"]

file= open("output.vtu","w+")
file.write("<?xml version=\"1.0\"?>\n")
file.write("<VTKFile type=\"UnstructuredGrid\">\n")
file.write("<UnstructuredGrid>\n")
file.write("<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n" %(numnodes,numel) )
file.write("<Points>\n")
file.write("<DataArray NumberOfComponents=\"3\" type=\"Float32\" format=\"ascii\" >\n")
for i in range (numnodes):
    file.write("%f\n %f 0.\n" %(node[i,0],node[i,1]))
file.write("</DataArray>\n")
file.write("</Points>\n")
file.write("<Cells>\n")
file.write("<DataArray Name=\"connectivity\" type=\"Int32\" format=\"ascii\" >\n")
for e in range (numel):
    for n in range (4):
        file.write("%d " %(elnodes.astype(int)[e][n]))
    file.write("\n")
file.write("</DataArray>\n")
file.write("<DataArray Name=\"offsets\" type=\"Int32\" format=\"ascii\" >\n")
offset=4
for e in range (numel):
    file.write("%d " %(offset))
    offset+=4
file.write("\n</DataArray>\n")
file.write("<DataArray Name=\"types\" type=\"UInt8\" format=\"ascii\" >\n")
for e in range (numel):
    file.write("10 ")
file.write("\n</DataArray>\n");
file.write("</Cells>\n")
file.write("<PointData Scalars=\"scalars\" format=\"ascii\">\n")

varinc=0    
for var in range(numvars): #Variables
    file.write("<DataArray Name=\" %s \" NumberOfComponents=\"%d\" type=\"Float32\" format=\"ascii\" >\n" %(varname[var],var_dim[var]))
    print ("var", var)
    for n in range (numnodes):
        print ("node ",n)
        imax=int(var_dim[var])
        for i in range(imax): 
            dof=int(varinc+var_dim[var]*n+i)
            print ("dof",dof)
            file.write("%f " %(Uglob[dof]))
            # print("vnrow",vnrow.astype(int))         
        file.write("\n")
    print ("varinc",varinc)
    varinc+=numnodes*var_dim[var]
    file.write("\n</DataArray>\n")

# v=0
# for n in range (numnodes):
    # for d in range (ndof):
        # print("v,Uglob[v]",v,Uglob[v])
        # file.write("%f " %(Uglob[v]))
        # v=v+1
    # file.write("\n")
# file.write("\n</DataArray>\n")
# file.write("<DataArray Name=\"Vel\" NumberOfComponents=\"2\" type=\"Float32\" format=\"ascii\" >\n")
# v=0
# for n in range (numnodes):
    # for d in range (2):
        # file.write("%f " %(vnxy[n,d]))
    # file.write("\n")
    # #i+=var_dim[0]
# file.write("\n</DataArray>\n")
file.write("</PointData>\n")
file.write("</Piece>\n")
file.write("</UnstructuredGrid>\n")
file.write("</VTKFile>\n")
file.close
