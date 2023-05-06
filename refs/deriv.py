from numpy import *
import numpy.matlib #for zeros

#This function only computes dEdUF and dEdUFvp, since the other are NULL
def calc_dEdU(Fd,Fvpd,NsigF,NFvp):
    #Fet_inv=linalg.inv(Fet)
    dEdU=[matrix(numpy.matlib.zeros((4, 16))),matrix(numpy.matlib.zeros((4, 20)))]
    #Eqn 1, comes from 
    #Fet[0][0]=0.5-0.5*(1.)
    TF=matrix([[0,0,0,1],[0,-1,0,0],[0,0,-1,0],[1,0,0,0]])
    dEdFed_inv=matrix(numpy.matlib.zeros((4, 5)))    #E.3
    
    #dFMdUF =arange(256).reshape(4,4,16)
    dFMdUF=numpy.zeros((4,4,16))

    Fet  = matrix(numpy.matlib.zeros((3,3)))#Tensor form    
    Fed=matrix(numpy.matlib.zeros((5, 1)))
    
    F4ed    =matrix(numpy.matlib.zeros((4, 1)))
    F4ed_inv=matrix(numpy.matlib.zeros((4, 1)))
    Fed_inv =matrix(numpy.matlib.zeros((5, 1)))
    
    F4vpd    =matrix(numpy.matlib.zeros((4, 1)))
    Fvpt     =matrix(numpy.matlib.zeros((3, 3)))
    Fvpt_inv =matrix(numpy.matlib.zeros((3, 3)))
    Fvpd_inv =matrix(numpy.matlib.zeros((4, 1)))
    
    F4vpd_inv=matrix(numpy.matlib.zeros((4, 1)))
    
    #Used for both F (E.7/E.8) 
    dFeinv_dFe=matrix(numpy.matlib.zeros((5, 5)))
    dF4einv_dF4e=matrix(numpy.matlib.zeros((4, 4)))
    
    dF4vpinv_dF4vp=matrix(numpy.matlib.zeros((4, 4)))
    
    dFedUF =matrix(numpy.matlib.zeros((5, 16)))
    dF4edUF=matrix(numpy.matlib.zeros((4, 16)))

    dFedUFvp =matrix(numpy.matlib.zeros((5, 20)))
    dF4edUFvp=matrix(numpy.matlib.zeros((4, 20)))

    ddetFe_dF4ed=matrix(numpy.matlib.zeros((1, 4)))
    ddetFvp_dF4vp=matrix(numpy.matlib.zeros((1, 4)))
    
    temp4=matrix(numpy.matlib.zeros((4, 1)))
    
    dF4vp_dUFvp     =matrix(numpy.matlib.zeros((4, 20)))
    
    #Check with inv
    #F=[Fxx]
    #print("Fd0",Fd[0,0])
    #E.11
    FM=matrix([[Fd[0,0],Fd[2,0],0    ,0   ],
               [Fd[1,0],Fd[3,0],0    ,0   ],
               [0    ,    0,Fd[0,0],Fd[2,0]],
               [0    ,    0,Fd[1,0],Fd[3,0]]])     
    
    #print ("FM",FM)
    
    # Ft[0,0]=Fd[0]
    # Ft[1,0]=Fd[1]
    # Ft[0,1]=Fd[2]
    # Ft[1,1]=Fd[3]
    # Ft[2,2]=1.

    for i in range(4):
        F4vpd[i]=Fvpd[i]
    
    Fvpt[0,0]=Fvpd[0]
    Fvpt[1,0]=Fvpd[1]
    Fvpt[0,1]=Fvpd[2]
    Fvpt[1,1]=Fvpd[3] 
    Fvpt[2,2]=Fvpd[4] 
    #print ("Fvpt",Fvpt)
    
    Fvpt_inv=linalg.inv(Fvpt) 
    F4vpd_inv[0]=Fvpt_inv[0,0]
    F4vpd_inv[1]=Fvpt_inv[1,0]
    F4vpd_inv[2]=Fvpt_inv[0,1]
    F4vpd_inv[3]=Fvpt_inv[1,1]
    #Fet=Ft*Fvpt_inv #Remains thermal part
    #print(Fet)
   
    #E.10 
    #TODO FOR FUTURE AND 3D FORMULATIONS:
    #Fe in fact is Fe=F*Fvp-1
    F4ed=FM*F4vpd_inv  
    #print ("Fvpd_inv",Fvpd_inv) 
    #Earr_e=
    #Earr_e=TLa*
    #E.9 
    for i in range (4):
        Fed[i]=F4ed[i]    
    Fed[4,0]=Fvpt_inv[2,2]       

    #print ("Fed(3)",Fed[3]) 
    
    det_Fed=Fed[0]*Fed[3]-Fed[1]*Fed[2] #E.6
    
    #e.5
    F4ed_inv=float(-1./det_Fed)*TF*F4ed
    
    for i in range(4):
        Fed_inv[i]=F4ed_inv[i]
    
    
    #Remember E_=[Exx Eyy Ezz Exy]
    #E.3 Pag 189
    dEdFed_inv[0,0]=dEdFed_inv[3,2]=Fed_inv[0]
    dEdFed_inv[0,1]=dEdFed_inv[3,4]=Fed_inv[1]
    dEdFed_inv[1,2]=dEdFed_inv[3,0]=Fed_inv[2]
    dEdFed_inv[1,3]=dEdFed_inv[3,1]=Fed_inv[3]
    dEdFed_inv[2,4]=Fed_inv[4]

    #E.14
    for k in range(8):
        dFMdUF[0,0,k]=dFMdUF[2,2,k]=NsigF[0,k]
        dFMdUF[0,1,k]=dFMdUF[2,3,k]=NsigF[2,k] 
        dFMdUF[1,0,k]=dFMdUF[3,2,k]=NsigF[1,k]
        dFMdUF[1,1,k]=dFMdUF[3,3,k]=NsigF[3,k] 
        
    #Third Term dFedUF
    #E.12
    for i in range(4):
        for j in range(16):
            for k in range (4):
                dF4edUF[i,j]=dFMdUF[i,k,j]*F4vpd_inv[k]
    
    #E.13
    for j in range(16):
        for i in range(4):
            dFedUF[i,j]=dF4edUF[i,j]
        dFedUF[4,j]=0.
    
    #Determinant Derivative, similar to C.22
    ddetFe_dF4ed[0,0]= F4ed[3,0]
    ddetFe_dF4ed[0,3]=-F4ed[2,0]
    ddetFe_dF4ed[0,1]=-F4ed[1,0]
    ddetFe_dF4ed[0,2]= F4ed[0,0]
    
    #E.7
    for i,k in range(4,4):
        temp4[i,0]=temp4[i,0]+TF[i,k]*F4ed[k,0]

    #print ("Fed",Fed)        
    #print ("det_Fed",det_Fed)
        
    dF4einv_dF4e=float(-1./(det_Fed*det_Fed))*temp4*ddetFe_dF4ed+float(1./(det_Fed))*TF
    
    #E.8
    for i,k in zip(range(4),range(4)):
        #print("i,j",i,j)
        dFeinv_dFe[i,k]=dF4einv_dF4e[i,k]
    
    dFeinv_dFe[4,4]=float(-1./(Fed[4,0]*Fed[4,0]))
        
    #print ("dF4einv_dF4e",dF4einv_dF4e)
    #F derivative
    #Eqn E.2
    dEdU[0]=dEdFed_inv*dFeinv_dFe*dFedUF
    
    #-----------------------------------------------------------------------------
    #Viscoplastic Fvp derivative
    #-----------------------------------------------------------------------------
    det_Fvp=Fvpd[0]*Fvpd[3]-Fvpd[1]*Fvpd[2]
    #E.17
    for i in range(4):
        temp4[i,0]=0.

    for i, k in zip(range(4), range(4)):
        temp4[i,0]=temp4[i,0]+TF[i,k]*F4vpd[k,0]

    ddetFvp_dF4vp[0,0]= F4vpd[3,0]
    ddetFvp_dF4vp[0,3]=-F4vpd[2,0]
    ddetFvp_dF4vp[0,1]=-F4vpd[1,0]
    ddetFvp_dF4vp[0,2]= F4vpd[0,0]
    
    dF4vpinv_dF4vp=float(-1./(det_Fvp*det_Fvp))*temp4*ddetFvp_dF4vp+float(1./(det_Fvp))*TF
    
    #First two terms are the same as before 
    #E.19
    for i, j in zip(range(4), range(20)):
        dF4vp_dUFvp[i,j]=NFvp[i,j]
    
    dF4edUFvp=FM*dF4vpinv_dF4vp*dF4vp_dUFvp
    
    #E.15
    for j in range(20):
        for i in range(4):
            dFedUFvp[i,j]=dF4edUFvp[i,j]
        dFedUFvp[4,j]=-1./(Fvpd[4]*Fvpd[4])*NFvp[4,j]
    
    #Eqn E.2
    dEdU[1]=dEdFed_inv*dFeinv_dFe*dFedUFvp
    
    #print("dEdU[0]",dEdU[0])
    #print("dEdU[1]",dEdU[1])
       
    return dEdU
