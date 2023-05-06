# https://docs.sympy.org/latest/tutorial/matrices.html
#Plain strain/stress 1 element example
#import math
from numpy import *
import numpy.matlib #for zeros

import array as arr
import deriv

#******************************************************
#----------------------------
#Input Data------------------
r0=1.
r1=5.
lx=r1-r0
ly=20.*3.1415926/180.
nex=20
ney=8
numit=20
solver=1 #1:simple 2:Newton Raphson

#RADIAL FLOW EXAMPLE 


#-------------



numvars=1 #1: Only F tensor, 2: F and internal variable s
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
mat_s0 =125.1e6

ey=1.e9
nu=0.499
##********************** INPUT END 

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

#BEFORE CONVERT TO CARTESIAN!
#Only for this example   
print("************VELOCITIES***************")
for n in range(numnodes):
    r=node[n,0]+r0
    t=node[n,1]-ly/2.
    vr=0.1*r0/r
    vnxy[n,0]=vr*cos(t)
    vnxy[n,1]=vr*sin(t)
    #print("vxy ",n,":",vnxy[n,0],vnxy[n,1])


#Radial flow example: convert from cilindrical to cartesian
for n in range(numnodes):
   r=node[n,0]+r0
   t=node[n,1]-ly/2.
   node[n,0]=r*cos(t)
   node[n,1]=r*sin(t)
   #print("Coord ",n,":",node[n,0],node[n,1])
 

#BEFORE BOUNDARY CONDITIONS
#Element dimension and DOF PER VARIABLE! 
var_edof=zeros(2)
var_dim =[4,1]
ndof=0
#Formulation type and DOFs
for i in range(numvars): 
    ndof += var_dim[i]
print ("dofs: ", ndof)
for i in range(numvars):
    var_edof[i]=4*var_dim[i]  
    
#### BOUNDARY CONDITIONS
# BOUNDARY CONDITIONS
boundarynode=zeros(ney+1)
#ROWS ARE NODES; COLS ARE ENTIRE BCs vector (all vars)
node_bc=matrix(numpy.matlib.zeros((size(boundarynode), ndof)))

dnode=(nex+1)    
i=0
for dy in range(ney+1): 
    inode=dy*dnode
    boundarynode[i]=inode
    print("i",i)
    node_bc[i]=[1.,0.,0.,1.]
    i+=1

print("node_bc",node_bc)
    
print("boundarynode",boundarynode)

################################################
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
print ("GLOBAL DOF: ",dof)
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
#BsigF=[matrix(numpy.matlib.zeros((4, 16))),matrix(numpy.matlib.zeros((4, 8)))]
BsigF=numpy.zeros((4,16,2))
BFvp =arange(200).reshape(5,20,2) #

temp4x16=matrix(numpy.matlib.zeros((4, 16)))

#(5,20,2) x (20)
BUFvp=matrix(numpy.matlib.zeros((5, 2)))

B4i=numpy.zeros((4,4,2))
B5i=numpy.zeros((5,5,2))

#(4,16,2)
#print(BsigF[0])
LM =matrix(numpy.matlib.zeros((4, 4)))
LM5=matrix(numpy.matlib.zeros((5, 5)))


K=matrix(numpy.matlib.zeros((44, 44)))

R   =[  matrix(numpy.matlib.zeros((16, 1))),
        matrix(numpy.matlib.zeros(( 4, 1)))]
Rzero =[  matrix(numpy.matlib.zeros((16, 1))),
         matrix(numpy.matlib.zeros(( 4, 1)))]
RF  =matrix(numpy.matlib.zeros((16, 1)))
Rsig=matrix(numpy.matlib.zeros((16, 1)))
RFvp=matrix(numpy.matlib.zeros((20, 1)))
Rs  =matrix(numpy.matlib.zeros((4, 1)))
Rv  =matrix(numpy.matlib.zeros((8, 1)))

#Formulation 1
#-----------------------------------
#Symmetric tensors
Ee =matrix(numpy.matlib.zeros((4, 1)))
Eet=matrix(numpy.matlib.zeros((2, 2))) #Tensor form 



#These are the same but reorganized
dVxy=zeros(4)
L   =matrix(numpy.matlib.zeros((2, 2)))
BL      = arange(128).reshape(4,4,8)            #Eqns 2.33, B.17
temp8x1 = matrix(numpy.matlib.zeros((8, 1))) 


#Stress
sig  =matrix(numpy.matlib.zeros((4, 1)))    #Stress Gauss Points
sig_d=matrix(numpy.matlib.zeros((4, 1)))    #Deviatoric, is also symmetric
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
 
#Matrix Double Entry
#Example in formulation 1:
# ( 8x 8)  (8x16)  (8x20) (8x4)
# (16x 8) (16 x 16) ..  ..
# (20x 8) (20x16) (20x20) (20x4)
# ( 4x 8)

dDdUv   =arange(200).reshape(5,5, 8) #TO MODIFY
dDdUF   =arange(400).reshape(5,5,16) #TO MODIFY
dDdUFvp =arange(500).reshape(5,5,20) #TO MODIFY
dDdUs   =arange(100).reshape(5,5, 4) #TO MODIFY

#From multiplying dDdU x NFvp*
#4.39 to 4.42
temp_dDFvp=[numpy.matlib.zeros((5, 8)),
            numpy.matlib.zeros((5,16)),
            numpy.matlib.zeros((5,20)),
            numpy.matlib.zeros((5, 4))]

Kt=[
     [matrix(numpy.matlib.zeros(( var_edof.astype(int)[0], var_edof.astype(int)[0]))), matrix(numpy.matlib.zeros((var_edof.astype(int)[0], var_edof.astype(int)[1])))]
     ,
     [matrix(numpy.matlib.zeros(( var_edof.astype(int)[1], var_edof.astype(int)[0]))), matrix(numpy.matlib.zeros((var_edof.astype(int)[1], var_edof.astype(int)[1])))]
    ] 
Kzero=[
     [matrix(numpy.matlib.zeros(( var_edof.astype(int)[0], var_edof.astype(int)[0]))), matrix(numpy.matlib.zeros((var_edof.astype(int)[0], var_edof.astype(int)[1])))]
     ,
     [matrix(numpy.matlib.zeros(( var_edof.astype(int)[1], var_edof.astype(int)[0]))), matrix(numpy.matlib.zeros((var_edof.astype(int)[1], var_edof.astype(int)[1])))]
    ] 
print ("Kt", Kt)
    
dgdU=[  matrix(numpy.matlib.zeros((1, 8))),
        matrix(numpy.matlib.zeros((1,16))),
        matrix(numpy.matlib.zeros((1,20))),
        matrix(numpy.matlib.zeros((1, 4)))]

print("Kt_0",len((Kt[1][0])),len((Kt[1][0]).transpose()))

ck = ey/ ((1. + nu)*(1. - 2. * nu))
c=matrix(numpy.matlib.zeros((4, 4)))
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

for n in range(numnodes):
    #Velocities 
    iF=ndof*n
    #Initial deformation gradients as identity??
    Uglob[iF  ]=Uglob[iF+3]=1
    Uglob[iF+1]=Uglob[iF+2]=0#xy and yx        

                        
print ("Initial Uglob", Uglob)

#-------------------------------------------
it=0

## ------------------------------------------
## Newton Rhapson Loop
## ------------------------------------------
while (it < numit):

    #Clean Global Matrices for assembly
    #print ("Kglob",Kglob)
    for idof in range(dof):
        Rglob [idof] = 0.
        for jdof in range(dof):
            Kglob[idof,jdof] = 0.    
#ELEMENT LOOP  ----------------
    #for e in range (1):
    for e in range (numel): # TO MODIFY 
        #Obtain Ve from global
        Kt=Kzero
        R=Rzero
        #print("Kt[0][0]",Kt[0][0])    
        #Obtain Ve from global
        Kel=0.
        for n in range(4):
            X2[n]=node[elnodes.astype(int)[e][n]]
            
        #print ("Element ", e)
        #print ("Element Nodes")
        #print ("X2",X2)
        #dHrs=matrix([[(1+sg),-(1+sg),-(1-sg),(1-sg)], [(1+rg),(1-rg),-(1-rg),-(1+rg)] ])
        dHrs=matrix([[(1+0.),-(1+0.),-(1-0.),(1-0.)], [(1+0.),(1-0.),-(1-0.),-(1+0.)] ])
        J=dHrs*X2
        lx=J[0,0]/2.
        ly=J[1,1]/2.
        #print("J,lx,ly",J,lx,ly)
                
        for ig in range(2):
            for jg in range(2):
                rg=gauss[ig]
                sg=gauss[jg]

                #Numerated as in Bathe
                Ns  =0.25*matrix([(1+sg)*(1+rg),(1-rg)*(1+sg),(1-sg)*(1-rg),(1-sg)*(1+rg)])   
                dHrs=matrix([[(1+sg),-(1+sg),-(1-sg),(1-sg)], [(1+rg),(1-rg),-(1-rg),-(1+rg)] ])
                #print ("Ns",Ns)
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
                
                #print("NsigF",NsigF)
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
      
                #print ("BsigF",BsigF)                
                #Interpolate velocity
                #INCREMENT GLOBAL VELOCITY FROM INCREMENTS!!!
                #CHANGE F TO ASSEMBLE IN THE SAME PLACE FOR BOTH FORMS
                juf=0
                iv=0
                for n in range (4): #Element nodes
                    gn=elnodes.astype(int)[e][n]#globnod
                    
                    for i in range (2):    #Velocity is var 0
                        #print("UV,i+iv,node",i+iv,gn)
                        UV[i+iv,0]=vnxy[gn,i]
                        
                    for j in range (var_dim[0]):
                        UF  [j+juf,0]=Uglob[ndof*gn+j]
                        #print("UF(j+juf,dof)",j+juf,ndof*gn+j)
                    juf+=var_dim[0]
                    iv=iv+2

                #print("UV",UV)
                #print("UF",UF)
                
                v  =Nv*UV #[2x8 x (8x1)]
                s  =float(Ns*Us)
                F  =NsigF*UF #[(4x16)*(16x1) =(4x1)]
                
                 
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
                
                #Stabilization factor tau 2.26
                #tau=beta*he/(2|v|)
                #See beta estability paramter
                #LM (2.33 & 4.23 p23 & p91)
                #Attention: SE DIFFERENCES WITH L_ in 2.28
                LM[0,0]=LM[2,2]=dVxy[0]
                LM[0,1]=LM[2,3]=dVxy[1]
                LM[1,0]=LM[3,2]=dVxy[2]
                LM[1,1]=LM[3,3]=dVxy[3]
                
                for i in range(4):
                    for j in range(4):               
                        LM5[i,j]=LM[i,j]
                    
                w=1. #TO MODIFY
                #Calculate sigma
                #2.31 Comes first from 2.2 (Then 2.31)
                #Pij=vk d(sig ij)/xk - dvi/dxk sig(kj) + dvk/dxk sig ij
                #Calculate Piola Kirchoff Pi (2.31) Gij Cjk-Gij LM (jk) sig(k)+
                #Attention double contraction
                #P=G*c*E-G*LM*sig+(LM[0,0]+LM[1,1])*G*sig
                beta=1.
                #print("u2+v2",v[0]*v[0]+v[1]*[1])
                vnorm=sqrt(v[0]*v[0]+v[1]*v[1])
                #he=(lx+ly)/2.
                he=(lx*v[0]+ly*v[1])/vnorm #HUGHES (APPROX)
                tau=float(beta*he/(2.*vnorm))
                #tau=0.
                #print ("tau",tau)
                #Calculate stabilization parameter
                #tau=1.
                #STRESSES**********
                #From 2.27 Plane Strain Symmetric tensors are defined as 
                #t=[txx tyy tzz tyz]
                pi=1./3.*(sig[0,0]+sig[1,0]+sig[2,0])
                for i in range(3): #Only daigonal is modified
                    sig_d[i,0]=sig[i,0]-pi #comps are [x y z yz]
                #print ("sigd",sig_d[i][0])   
                for k in range(4):
                    sig_eq=sqrt(1.5*(sig_d[k,0]))
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

                #for i in range(2):
                #    for j in range(2):
                #Dvp comes from N, which comes from 
                #Dvp is [x y z yz]
                DM[0,0]=DM[2,2]=Dvp[0]
                DM[1,1]=DM[3,3]=Dvp[1]
                DM[1,0]=DM[3,2]=DM[0,1]=DM[2,3]=Dvp[3]
                DM[3,3]=Dvp[2]
                
                wJ=w*detJ
                
             
                visc=1.
       
                # *****************************************************************
                #RESIDUALS ******************* 2.26 to 2.39 *****************************
                #Construct vk Bsig mik
                for m in range(4):
                    for i in range(16):
                            temp4x16[m,i]=0.
                for m in range(4):
                    for i in range(16):
                        for k in range(2):
                            temp4x16[m,i]=temp4x16[m,i]+BsigF[m,i,k]*v[k,0]
               
                
                R[0]   = R[0] + (NsigF+temp4x16*tau).transpose()*(temp4x16*UF-LM*NsigF*UF)*wJ
                if numvars == 2:
                    R[1]    = R[1] + (Ns+tau*v.transpose()*Bs).transpose()*(v.transpose()*Bs*Us-g_sigs)*wJ
              
                #R Assembly            
                #TANGENT MATRIX   
                #PAGES 25 y 94
                
                #----------------------------------------------------                   
                               
                #dRFdUF  4.36

                Kt[0][0]= Kt[0][0]+(
                                (NsigF+temp4x16*tau).transpose()*
                                (temp4x16-LM*NsigF)
                                )*wJ
                               
                #print("temp4x16-LM*NsigF",temp4x16-LM*NsigF)         
                #print("temp4x16",temp4x16)                          
                #dRFdUF=dRFdUF=0.  4.37 & 4.38                
                            
                #---------------------------------------------------
                #S derivatives -- COMMON TO BOTH FORMULATIONS
                #2.53 and 4.xx ATENTION IN CHAPTER 4 k and p indices are wrong, see 2.53
                #FORMER 3,1 and 3,3

                #4.44 dRs/dF 4.44 (4x16)
                if numvars == 2 :
                    #print("Kt(1,0)",Kt[1][0])
                    Kt[1][0]=Kt[1][0]+(
                              (Ns+tau*v.transpose()*Bs).transpose()*
                              (-dgdU[1])*
                              wJ) 
                    
                    #dRs/dS 4.46    
                    Kt[1][1]=Kt[1][1]+(
                              (Ns+tau*v.transpose()*Bs).transpose()*
                              (v.transpose()*Bs-dgdU[3])*
                              wJ) 
        
        #print ("Nv",Nv)
        #print("Kt[0][0]",Kt[0][0]) 
        
        vrowinc=0
        #Assembly Matrix
        for vrow in range(numvars): #Variables
            ir=0
            imax=int(var_dim[vrow])
            for n in range (4): #Nodes
                for i in range(imax): 
                    d=elnodes.astype(int)[e][n]
                    #print("vrowinc,d,a+b",vrowinc,d,vrowinc+var_dim[vrow]*d+i)
                    vnrow[ir]=vrowinc+var_dim[vrow]*d+i
                    ir=ir+1
            
            #print ("vrow vnrow",vrow,vnrow)
            vcolinc=0        
            for vcol in range(numvars): #Variables
                
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
                
                #print("vcol vncol",vcol,vncol)
                for row in range(4*imax):
                    if solver==2:   #NEWTON RAPHSON
                        Rglob[vnrow.astype(int)[row]]-=R[vrow][row]
                    for col in range(4*jmax):
                        #print("vnrow(row)vncol(col)",vnrow[row],vncol[col]) 
                        Kglob[vnrow.astype(int)[row],vncol.astype(int)[col]] =  Kglob[vnrow.astype(int)[row],vncol.astype(int)[col]]+(
                                                                              Kt[vrow][vcol][row,col])
                vcolinc+=numnodes*var_dim[vcol]
            
            
            
            vrowinc+=numnodes*var_dim[vrow]
        

    ##---------------------------------------------------------------
    ##Boundary conditions
    ##---------------------------------------------------------------
    #Velocity DOFs 
    #In this example velocities are known
    #AT INLET(left nodes):
    # F=I , sigma = 0
    if solver == 1: #NONZERO VALUES
        for n in range(size(boundarynode)):
            for i in range ( var_dim [ 0 ] ):
                inode=boundarynode[n]
                idof = var_dim[0] * inode + i
                print("idof",idof)
                for j in range(dof):
                    Rglob[ j ] = Rglob[ j ] - Kglob[j,int(idof)] * node_bc[ n, i ] #dU=0, U=1(idof)
            
    dnode=(nex+1)    
    for dy in range(ney+1): 
        inode=dy*dnode
        #print("node",inode)   
        #Deformation gradient F
        for i in range ( var_dim [ 0 ] ):
            idof = var_dim[0] * inode + i
            #print ("idof",idof)
            for j in range(dof):
                Kglob[ idof , j ] = 0
                #Rglob[ j ] -= Kglob[j,idof] * 0 #dU=0, U=1(idof)
                Kglob[ j ,idof ] = 0
            

            Kglob[idof,idof] = 1
            if solver == 2:
                Rglob[idof  ] = 0           #F INCREMENT (dF) IS NULL!!!!!
        
        #Sigma is zero, Internal variable s ,      
        if numvars == 2:
            idofs = idof + var_dim[0]
            for i in range(dof):                    
                Kglob[ idofs , i ] = 0
                #Rglob[ i ] -= Kglob[ i, idofs ] * 0 #1 is R(idof)  
                Kglob[ i , idofs ] = 0                
        
            Kglob[idofs,idofs] = 1
            Rglob[idofs  ] = 0                
    
    #BOUNDARY CONDITIONS IN STANDARD SOLVER (SUCCESIVE ITERATIONS)
    if solver == 1:
        for n in range(size(boundarynode)):
            inode=boundarynode[n]
            for i in range ( var_dim [ 0 ] ):
                idof = var_dim[0] * inode + i   
                print("i, idof",i, idof)
                Rglob[ int(idof) ] = node_bc[ n, i ] #dU=0, U=1(idof)      
                
    # print("KGLOB\n")
    # for i in range (dof):
        # for j in range (dof):
            # print(Kglob[i,j], end = " ")
        # print("\n")   
        # # print("\n")
    # print("Rglob",Rglob)


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
    print ("Iteration: ",it, "from ", numit)
    


    
#print ("Results")
#print("Uglob", Uglob)

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
file.write("<DataArray Name=\"DefGrad_F\" NumberOfComponents=\"%d\" type=\"Float32\" format=\"ascii\" >\n" %(ndof))
v=0
for n in range (numnodes):
    for d in range (ndof):
        print("v,Uglob[v]",v,Uglob[v])
        file.write("%f " %(Uglob[v]))
        v=v+1
    file.write("\n")
file.write("\n</DataArray>\n")
file.write("<DataArray Name=\"Vel\" NumberOfComponents=\"2\" type=\"Float32\" format=\"ascii\" >\n")
v=0
for n in range (numnodes):
    for d in range (2):
        file.write("%f " %(vnxy[n,d]))
    file.write("\n")
    #i+=var_dim[0]
file.write("\n</DataArray>\n")
file.write("</PointData>\n")
file.write("</Piece>\n")
file.write("</UnstructuredGrid>\n")
file.write("</VTKFile>\n")
file.close
