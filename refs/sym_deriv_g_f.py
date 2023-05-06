from sympy import *
#from numpy import *
import numpy.matlib #for zeros
import numpy as np

#From 4.2 
#sig= C Ee + 2 mu D
#From 4.3 and 4.4 
#Ee = 1/2 (deltaij -Fe(-1)Fe(-1) ) 4.3 
#Dij= 1/2 (Lji + Lij )
#From 4.6 
#Fe=F Ftem-1 Fvp-1


s,h0,A,s_,s_ast,m,n,psi,Uf,Ufvp,Dvp,a = symbols('s h0 A s_ s_ast m n psi Uf Ufvp Dvp a')

#From evp=f(s,sig)

z = symbols('z', cls=Function)
Ee = symbols('Ee', cls=Function)
sig=z(Uf,Ufvp)
f= A*Pow(psi*sinh(s/sig),1/m)
param=1 - s/(s_*Pow(psi*sinh(s/sig),1/m))

#Cannot anidate 2 pows and differenciate
#(1-s/s)a*=(1-s/(s_*( A (sinh)  /A)n)
#1-s/s* > 0
g = h0*(Pow(param,a))*Pow(f,1/m)


#--------- DVP
sig_t=np.matrix(np.matlib.zeros((4, 1)))    #Deviatoric, is also symmetric
sig_d=np.matrix(np.matlib.zeros((4, 1)))    #Deviatoric, is also symmetric

#pi=1./3.*(sig[0,0]+sig[1,0]+sig[2,0])
#for i in range(3): #Only daigonal is modified
#    sig_d[i,0]=sig[i,0]-pi #comps are [x y z yz]
#print ("sigd",sig_d[i][0])   
#for k in range(4):
#    sig_eq=sqrt(1.5*(sig_d[k,0]))
#Eq 4.9

#Dvp=matrix(numpy.matlib.zeros((4,4)))
#sqrt(3/2)*f


#1-s/s* < 0
g2 = h0*(1-s/(s_*Pow(f/A,n)))*f

gt=Pow(1-s,a)
print(simplify(diff(gt,s)))

print("f's: ",diff(f,s))
print("f'Uf: ",diff(f,Uf))
print("-------------------------")
print("**********************************************")
print("(1-s/s*)>0")
print("g1's NOT SIMPL : ",diff(g,s))
print("******************************")
print("g1's SIMPLIFIED: ",simplify(diff(g,s)))
print("-------------------------")
print("g1'Uf NOT SIMP: ",diff(g,Uf))
print("******************************")
print("g1'Uf NOT SIMP: ",simplify(diff(g,Uf)))
print("-------------------------")

from sympy import *
r, t = symbols('r t') # r (radius), t (angle theta)
f = symbols('f', cls=Function)
x = r* cos(t)
y = r* sin(t)
g = f(x,y)
#print(Derivative(g,r, 2).doit())

