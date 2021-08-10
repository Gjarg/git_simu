import numpy as np


class Crystal:
    def __init__(self, name, theta, phi=None, plan=None):
        self.name = str(name)
        self.theta = theta
        self.phi = phi
        self.plan = plan