from d2q9 import D2Q9


class BoundaryConditions(D2Q9):

    def __init__(self, f):
        self.f = f
        self.f[:, 0, :] = self.f[:, 1, :]
        self.f[:, -1, :] = self.f[:, -2, :]
        self.f[:, :, 0] = self.f[:, :, 1]
        self.f[:, :, -1] = self.f[:, :, -2]
        D2Q9.__init__(self)
        
    def periodic_in_x(self):
        """
        To apply periodic boundary condition in x
        """
        self.f[1, 0, 1:-1] = self.f[1, -1, 1:-1]
        self.f[8, 0, 1:-1] = self.f[8, -1, 1:-1]
        self.f[5, 0, 1:-1] = self.f[5, -1, 1:-1]

        self.f[6, -1, 1:-1] = self.f[6, 0, 1:-1]
        self.f[3, -1, 1:-1] = self.f[3, 0, 1:-1]
        self.f[7, -1, 1:-1] = self.f[7, 0, 1:-1]

    def bounce_back_right(self):
        """
        Apply bounce back on right wall
        """
        self.f[3, -2, 2:-2] = self.f[1, -1, 2:-2]
        self.f[6, -2, 2:-2] = self.f[8, -1, 2:-2]
        self.f[7, -2, 2:-2] = self.f[5, -1, 2:-2]

    def bounce_back_left(self):
        """
        Apply bounce back on right wall
        """
        self.f[8, 1, 2:-2] = self.f[6, 0, 2:-2]
        self.f[1, 1, 2:-2] = self.f[3, 0, 2:-2]
        self.f[5, 1, 2:-2] = self.f[7, 0, 2:-2]

    def bounce_back_bottom(self):
        """
        Apply bounce back on right wall
        """
        self.f[2, 2:-2, 1] = self.f[4, 2:-2, 0]
        self.f[6, 2:-2, 1] = self.f[8, 2:-2, 0]
        self.f[5, 2:-2, 1] = self.f[7, 2:-2, 0]

    def moving_wall_bc_top(self, u_wall):
        """
        Apply moving wall bounce back on top wall
        """
        density=self.f.sum(axis=0)
        self.f[4, 1:-1, -2] = self.f[2, 1:-1, -1]
        self.f[7, 1:-1, -2] = self.f[5, 1:-1, -1] + 6.0 * density[1:-1, -1] * self.w[7] * self.c[0][7] * u_wall
        self.f[8, 1:-1, -2] = self.f[6, 1:-1, -1] + 6.0 * density[1:-1, -1] * self.w[8] * self.c[0][8] * u_wall

    def bot_left_corner_correction(self):
        """
        Bottom left corner
        """
        self.f[2, 1, 1] = self.f[4, 1, 0]
        self.f[5, 1, 1] = self.f[7, 1, 0]
        self.f[1, 1, 1] = self.f[3, 1, 0]
        self.f[8, 1, 1] = self.f[6, 1, 0]
        self.f[6, 1, 1] = self.f[8, 1, 0]

    def bot_right_corner_correction(self):
        """
        Bottom right corner
        """
        self.f[2, -2, 1] = self.f[4, -2, 0]
        self.f[6, -2, 1] = self.f[8, -2, 0]
        self.f[3, -2, 1] = self.f[1, -2, 0]
        self.f[5, -2, 1] = self.f[7, -2, 0]
        self.f[7, -2, 1] = self.f[5, -2, 0]

    def top_left_corner_correction(self):
        """
        Top left corner
        """
        self.f[5, 1, -2] = self.f[7, 1, -1]
        self.f[7, 1, -2] = self.f[5, 1, -1]
        self.f[1, 1, -2] = self.f[3, 1, -1]
        self.f[8, 1, -2] = self.f[6, 1, -1]
        self.f[4, 1, -2] = self.f[2, 1, -1]