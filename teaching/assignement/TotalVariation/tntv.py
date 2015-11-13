# Initialization
from dolfin import *
import math
import numpy as np
import logging
import matplotlib.pyplot as plt
import nb
from unconstrainedMinimization import InexactNewtonCG

logging.getLogger('FFC').setLevel(logging.WARNING)
logging.getLogger('UFL').setLevel(logging.WARNING)
set_log_active(False)


# Set the level of noise:
noise_std_dev = .3

# Load the image from file
data = np.loadtxt('image.dat', delimiter=',')
np.random.seed(seed=1)
noise = noise_std_dev*np.random.randn(data.shape[0], data.shape[1])

# Set up the domain and the finite element space.
Lx = float(data.shape[1])/float(data.shape[0])
Ly = 1.

mesh = RectangleMesh(0,0,Lx,Ly,200, 100)
V = FunctionSpace(mesh, "Lagrange",1)

# Generate the true image (u_true) and the noisy data (u_0)
class Image(Expression):
    def __init__(self, Lx, Ly, data):
        self.data = data
        self.hx = Lx/float(data.shape[1]-1)
        self.hy = Ly/float(data.shape[0]-1)
        
    def eval(self, values, x):
        j = math.floor(x[0]/self.hx)
        i = math.floor(x[1]/self.hy)
        values[0] = self.data[i,j]

trueImage = Image(Lx,Ly,data)
noisyImage = Image(Lx,Ly,data+noise)
u_true  = interpolate(trueImage, V)
u_0 = interpolate(noisyImage, V)

nb.multi1_plot([u_true,u_0], titles=["True Image", "Noisy Image"])
plt.show()

# 

