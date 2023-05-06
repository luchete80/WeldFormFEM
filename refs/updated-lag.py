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
nex=40
ney=5
plastic=0
plastic=0
numit=3
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



dVxy=zeros(4)
L   =matrix(numpy.matlib.zeros((2, 2)))
BL      = numpy.zeros((3,16))                      #BATHE, pp. 
BNL      = numpy.zeros((4,16))                      #BATHE, pp. 
temp8x1 = matrix(numpy.matlib.zeros((8, 1))) 


#Stress
sig  =matrix(numpy.matlib.zeros((4, 1)))    #Stress Gauss Points xx, yy,zz, syx, is sbar, WHICH IS SYMMETRIC
sig_d=matrix(numpy.matlib.zeros((4, 1)))    #Deviatoric, is also symmetric
signod=matrix(numpy.matlib.zeros((numnodes, 4)))
sigeq_nod=zeros(numnodes)
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
    for k in range(4): #Nodes
        #shape functions
        Nv[0,2*k  ]=Nv[1,2*k+1]=Ns[0,k]
        for j in range(4):
            NsigF[j,4*k+j]=Ns[0,k]

        #derivatives Bv (B.14)
        Bv[0,2*k  ]=dHxy[0,k]
        Bv[1,2*k  ]=dHxy[1,k]
        Bv[2,2*k+1]=dHxy[0,k]
        Bv[3,2*k+1]=dHxy[1,k]
        
        BL[0,2*k  ]=dHxy[0,k]
        BL[1,2*k+1]=dHxy[1,k]
        BL[2,2*k  ]=dHxy[1,k]
        BL[2,2*k+1]=dHxy[0,k]
        

        BNL[0,2*k  ]=dHxy[0,k]
        BNL[1,2*k  ]=dHxy[1,k]
        BNL[2,2*k+1]=dHxy[0,k]
        BNL[3,2*k+1]=dHxy[1,k]
        
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

    return Nv,NsigF,Ns,BL,BsigF,BFvp

def calculate_Jacobian(rg,sg,X2):
    dHrs=matrix([[(1+sg),-(1+sg),-(1-sg),(1-sg)], [(1+rg),(1-rg),-(1-rg),-(1+rg)] ])
    #Numerated as in deal.ii
    #dHrs=matrix([[-(1-s),(1-s),-(1+s),(1+s)], [-(1-r),-(1+r),(1-r),(1+r)] ])        
    dHrs/=4.
    J=dHrs*X2
    dHxy=linalg.inv(J)*dHrs
    detJ=linalg.det(J) 
    return J,detJ
    
def calculate_Uelem(numnodes,elnodes,Uglob,var_dim):
    if plastic ==0:
        numvars=2
    else:
        numvars=4
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
    if plastic ==0:
        return UV,UF
    else:
        return UV,UF,UFvp,Us
        
def calculate_Vderivs(Bv,UV):
    dVxy=Bv*UV #(4x8)*(8x1)=(4x1) (vx,x vx,y vy,x vy,y)T 
    L[0,0]=dVxy[0]
    L[0,1]=dVxy[1]
    L[1,0]=dVxy[2]
    L[1,1]=dVxy[3]
    
    Dt=0.5*(L+L.transpose())
    D[0]=Dt[0,0];D[1]=Dt[1,1];
    D[3]=Dt[0,1];
    return dVxy,L,Dt,D
    
def calculate_sigma(c,D):

    #print("F",F)
    #Ft=
    if (form==1):
        sig=NsigF*Usig
    else:
        if plastic:
            Fvp=NFvp*UFvp
    

    Ft[0,0]=F[0]
    Ft[0,1]=F[1]
    Ft[1,0]=F[2]
    Ft[1,1]=F[3]
    Ft[2,2]=1.      #Fzz, plain strain 
    
    
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
                    
    return sig,sig_eq




    
#CALCULATE NEIGHBOURS NODES
shared_nodes=zeros(numnodes)
for e in range (numel):
    for n in range(4):
        en=elnodes.astype(int)[e][n]
        shared_nodes[en]+=1

print ("shared_nodes",shared_nodes)

#STRESS CALCULATIONS (AVERAGED)
nod_gaussp=numpy.matrix([[1., 1.],[-1., 1.],[-1.,-1.],[1.,-1.]])
for e in range (numel):
    for n in range(4):
        X2[n]=node[elnodes.astype(int)[e][n]]
    
    print ("X2",X2)
    for i in range(4):
        rg=nod_gaussp[i,0];sg=nod_gaussp[i,1]
        print ("rg,sg",rg,sg)
        Nv,NsigF,Ns,BL,BsigF,BFvp=calculate_shape_matrices (rg, sg,X2)
        J,detJ=calculate_Jacobian(rg,sg,X2)
        if plastic==0:
            UV,UF=calculate_Uelem(numnodes,elnodes,Uglob,var_dim)
            
        dVxy,L,Dt,D=calculate_Vderivs(Bv,UV)
        print ("D",D)

        Ft[0,0]=F[0]
        Ft[0,1]=F[1]
        Ft[1,0]=F[2]
        Ft[1,1]=F[3]
        Ft[2,2]=1.      #Fzz, plain strain 
        
        
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
        print ("sig",sig)
        #signod[n]+=sig
        sigeq_nod[n]+=sig_eq
        #sigeq_nod
#average
for n in range (numnodes):
    #signod[n]/=shared_nodes[n]
    sigeq_nod[n]/=shared_nodes[n]
    

print("sig,sig_eq",signod,sigeq_nod)

# print ("nod_gauss",nod_gauss)

    
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
            #print ("dof",dof)
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
