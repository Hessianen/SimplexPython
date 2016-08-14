class Voigt(object):
    'A Voigt-profile curve object with a set of parameters and the calculated y-values from the parameters and the x-vector'

    def __init__(self,x,beta):
        self.nu  = beta[0]
        self.w_L = beta[1]
        self.w_G = beta[2]
        self.x_0 = beta[3]
        self.L_0 = beta[4]
        self.G_0 = beta[5]
        self.L = np.divide(L_0,(1+(np.divide(x-x_0,w_L))**2))
        self.G = G_0*np.exp(-m.log(2)*(np.divide(x-x_0,w_G))**2)
        self.V = nu*L+(1-nu)*G
