import numpy as np


class D2Q9:

    def __init__(self):
        self.c = np.array([[0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0],
                           [0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0]])
        self.w = np.array([4.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0,
                           1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0])
        self.q = 9
        self.T0 = 1 / 3.0
        self.a = np.sqrt(self.T0)

    def feq(self, moments):
        """
        To compute equilibrium distribution function.
        """
        
        feq = np.zeros((self.q, moments.shape[1], moments.shape[2]))
            
        udotc = (np.array([[self.c[0]]]).T*moments[1]) + (np.array([[self.c[1]]]).T*moments[2])
        usq = moments[1] * moments[1] + moments[2] * moments[2]
        feq = (np.array([[self.w]]).T*moments[0]) * (1.0 + 3.0 * udotc - 1.5 * usq + 4.5 * udotc ** 2
                                             + 4.5 * udotc ** 3 - 4.5 * usq * udotc)
        return feq

    def moments(self, f):
        """
        Compute moments (rho, u, v) from populations 
        """
        mom = np.zeros((3, f.shape[1], f.shape[2]))

        mom[0] = f.sum(axis=0)
        mom[1] = (np.array([[self.c[0]]]).T*f).sum(axis=0)/mom[0]
        mom[2] = (np.array([[self.c[1]]]).T*f).sum(axis=0)/mom[0]

        return mom

    def collision(self, f, beta, alpha=2.0):
        """
        computes collision term 
        """

        feq = self.feq(self.moments(f))
        
        f += alpha * beta * (feq - f)
        f.cumsum(axis=0)
        

    @staticmethod
    def advection(f):
        """
        Streaming of velocities
        """

        f[3, :-1, :] = f[3, 1:, :]
        f[4, :, :-1] = f[4, :, 1:]
        f[7, :-1, :-1] = f[7, 1:, 1:]
        f[8, 1:, :-1] = f[8, :-1, 1:]

        f[1, -1:0:-1, ::-1] = f[1, -2::-1, ::-1]
        f[2, ::-1, -1:0:-1] = f[2, ::-1, -2::-1]
        f[5, -1:0:-1, -1:0:-1] = f[5, -2::-1, -2::-1]
        f[6, -2:0:-1, -1:0:-1] = f[6, -1:1:-1, -2::-1]
