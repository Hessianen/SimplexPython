import math as m
import numpy as np
import matplotlib.pyplot as plt
import time as t

def voigt(x,beta):
    nu  = beta[0]
    w_L = beta[1]
    w_G = beta[2]
    x_0 = beta[3]
    L_0 = beta[4]
    G_0 = beta[5]  
    L = np.divide(L_0,(1+(np.divide(x-x_0,w_L))**2))
    G = G_0*np.exp(-m.log(2)*(np.divide(x-x_0,w_G))**2)
    V = nu*L+(1-nu)*G
    return V
    
def err(y,V):
    N = np.array(y-V)
    S = (np.linalg.norm(N))**2
    return S

s_pnt = 400                         #startpoint
e_pnt = 410                         #endpoint
n_pnt = np.divide(e_pnt-s_pnt,0.05) + 1      #number of points 
n_pnt = m.trunc(n_pnt)              #Convert from float to an integer
x = np.linspace(s_pnt,e_pnt,n_pnt)   #Create the x-axis with equal amount of points as the y-axis (since it is a IGOR-wave)

#Import the textfile to data
y = np.loadtxt("") #Path to data
y = np.array(y)
y = np.flipud(y)


offset = 250
y = y - offset
#line, = plt.plot(x, y, '-', linewidth=2)
#plt.show()

#The model part starts here. The profile will be a pseudo-Voigt (not a real
#convolution)

#We need three initial guesses in the simplex algorithm - but we start with
#one and then just shift the initial one with a positive/negative change in the
#parameter-vector.

#Three inital guesses of the parametervector beta for the pseudo-Voigt
#profile
nu = 0.2
w_L = 0.1
w_G = 0.5
x_0 = 405.7
L_0 = 700
G_0 = 700
#beta is the parameter tensor
beta = np.zeros((7,6));
beta[0,:] = [nu,w_L,w_G,x_0,L_0,G_0]
delta = [1.01,1.01,1.01,1.01,1.01,1.01,1.01]
for i in range(0, 6):
    beta[i+1,:] = beta[0,:]
    beta[i+1,i] = beta[0,i]*delta[i]
    
V = np.zeros((7,n_pnt))
S = np.zeros(7)
S = np.array(S)

for i in range(0,7):
        V[i,:] = voigt(x,beta[i,:])
        S[i] = err(y,V[i,:])

n = np.argsort(S)
S = np.sort(S)

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(x,y,'--', linewidth=2)
for i in range(0,7):
    ax1.plot(x,V[i,:],'-', linewidth=2)
plt.show()

for k  in range(0,100):
    t.sleep(0.1)
    plt.cla()
    ax1.plot(x,y,'--', linewidth=2)
    for i in range(0,7):
        ax1.plot(x,V[i,:],'-', linewidth=2)
    plt.show()
    for i in range(0, 7):
        V[i,:] = voigt(x,beta[i,:])
        S[i] = err(y,V[i,:])
    
    n = np.argsort(S)
    S = np.sort(S)

    
    #l,s,h represent the best, second worse and the worse case respectively
    S_l = S[0]
    S_s = S[5]
    S_h = S[6]
    
    beta_l = beta[n[0],:]
    beta_s = beta[n[5],:]
    beta_h = beta[n[6],:]
    
    #beta_m is the mean value of all betas except the worst, i.e. h
    beta_m = np.zeros(6)
    for i in range(0,6):
        beta_m = beta_m + beta[n[i],:]
    beta_m = np.divide(beta_m,6)
    
    beta_r = 2*beta_m - beta_h
    V_r = voigt(x,beta_r)
    S_r = err(y,V_r)
    beta_e = 2*beta_r - beta_m
    V_e = voigt(x,beta_e)
    S_e = err(y,V_e)

    if (S_l <= S_r) and (S_r < S_s):
        beta[n[6],:] = beta_r  #Reflection is the best
        #print('Reflect')
    elif (S_e < S_l):
        beta[n[6],:] = beta_e  #Expand is the best
        #print('Expand')
    elif (S_r >= S_s):
        if (S_r < S_h):
            beta_c = beta_m + 0.5*(beta_r-beta_m)   #Outside contraction is the better choise
            V_c = voigt(x,beta_c)
            S_c = err(y,V_c)
            if (S_c <= S_r):
                #print('Outside contraction')
                beta[n[6],:] = beta_c
            else:
                for i in range(1,7):
                    beta[i,:] = beta[n[0],:] + 0.5*(beta[n[i],:]-beta[n[0],:])
                #print('Shrink1')
        else:
            beta_c = beta_m + 0.5*(beta_h-beta_m)
            V_c = voigt(x,beta_c)
            S_c = err(y,V_c)
            if (S_c < S_h):
                #print('Inside contraction')
                beta[n[6],:] = beta_c
            else:
                for i in range(1,7):
                    beta[i,:] = beta[n[0],:] + 0.5*(beta[n[i],:]-beta[n[0],:])
                #print('Shrink2')
    else:
        for i in range(1,7):
            beta[i,:] = beta[n[0],:] + 0.5*(beta[n[i],:]-beta[n[0],:])
        #print('Shrink3')

#fig = plt.figure()
#plt.cla()
#ax1 = fig.add_subplot(111)
#ax1.plot(x,y,'--', linewidth=2)
#for i in range(0,7):
#    ax1.plot(x,V[i,:],'-', linewidth=2)
#plt.show()

print('nu = ')
print(beta[n[1],0])
print('w_L = ')
print(beta[n[1],1])
print('w_G = ')
print(beta[n[1],2])
print('x_0 = ')
print(beta[n[1],3])
print('L_0 = ')
print(beta[n[1],4])
print('G_0 = ')
print(beta[n[1],5])
