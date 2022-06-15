from d2q9 import D2Q9


class SimParams(D2Q9):
    def __init__(self, mach, reynolds_num, ref_len, res):
        self.mach_num = mach
        self.reynolds = reynolds_num
        self.ref_len = ref_len
        self.res = res
        self.dt = self.res
        D2Q9.__init__(self)

    @property
    def ref_vel(self):
        """
        Compute reference velocity
        """
        return self.mach_num * self.a

    @property
    def knudsen(self):
        """
        Compute Knudsen number
        """
        return self.mach_num / self.reynolds

    @property
    def nu(self):
        """
        Compute viscosity from reference length, velocity and Reynolds number
        """
        return self.ref_len * self.ref_vel / self.reynolds

    @property
    def tau0(self):
        return self.nu / self.T0

    @property
    def tau_ndim(self):
        return self.tau0 / self.dt

    @property
    def beta(self):
        """
        Relaxation 
        """
        return 1.0 / (2.0 * self.tau_ndim + 1.0)
