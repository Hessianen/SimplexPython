#err calculates the difference between the data y and the sum of the V's point by point, returns the norm of the difference to the power of 2.
def err(y,V):
    N = np.array(y)
    for i in V:
        N = N - np.array(V[i]) #Calculates the difference (reside) between the data in variable y and the peak V stored in an numpy array
    S = (np.linalg.norm(N))**2  #S is the vectornorm to power of 2
    return S
