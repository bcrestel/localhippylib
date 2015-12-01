import numpy as np
import socket
from time import sleep

from fenics import Vector
from posterior import *
from randomizedEigensolver import *
from variables import *
from NewtonCG import ReducedSpaceNewtonCG
from reducedHessian import ReducedHessian
from linalg import to_dense

from scipy.linalg import eigh

class MyOperator:
    def __init__(self, op,op2):
        self.op = op
        self.op2 = op2
        self.op2x = Vector()
        self.z = Vector()
        self.op.init_vector(self.op2x,1)
        self.op.init_vector(self.z,0)
            
    def init_vector(self,x,dim):
        self.op2.init_vector(x,1)
            
    def mult(self,x,y):
        self.op2.mult(x,self.op2x)
        self.op.solve(self.z,self.op2x)
        self.op2.mult(self.z,y)

class Server:
    def __init__(self, model):
        self.model = model
        TCP_IP = '127.0.0.1'
        TCP_PORT = 1983
        self.sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        self.sock.bind((TCP_IP, TCP_PORT))
        
        self.x_map = []
        self.d_misfit = []
        self.U_misfit = []
        self.d_gaussianPost = []
        self.U_gaussianPost = []
        self.laplace_approx_posterior = []
        self.inner_tol = 1e-15
        
        self.is_started = False
        
    def start(self):
        self.is_started = True
        self.sock.listen(1)
        
        print "waiting for connection"
        conn, addr = self.sock.accept() # hangs until other end connects
        print 'Connection address:', addr
        
        while True:
            print "waiting for message"
            cmd_class = conn.recv(32).strip()
            print "cmd_class: ",cmd_class
            if cmd_class == 'ComputeMapPoint':
                err = self.ComputeMapPoint()
                print "Sending", err
                if err:
                    conn.send("0")
                else:
                    conn.send("1")
                    
            elif cmd_class == 'KLE_GassianPost':
                k = self.KLE_GaussianPost()
                print "Sending", k
                conn.send("{0:d}".format(k))
            elif cmd_class == 'Eval':
                shape_str = conn.recv(30).strip()
                print "shape_str: ",shape_str
                shape = tuple(map(int,shape_str.split()))
                n_floats = np.prod(shape)
                eta_flat = np.zeros(n_floats,'float')        
                n_read = 0
                while n_read < n_floats:
                    n_toread = min(128,n_floats-n_read)
                    arr = np.fromstring(conn.recv(n_toread*8),dtype='>f8')
                    eta_flat[n_read:n_read+n_toread] = arr
                    n_read += n_toread
                    print "%i/%i floats read"%(n_read,n_floats)           
                eta = eta_flat.reshape(shape,order="F")
                print "received array: ",eta
                val = self.Eval(eta)
                print "send array: ",val
                conn.send("{0:1.16e}".format(val))
            elif cmd_class == 'quit':
                self.quit()
                break
            
    def ComputeMapPoint(self):
        a0 = self.model.prior.mean.copy()
        solver = ReducedSpaceNewtonCG(self.model)
        solver.parameters["rel_tolerance"] = 1e-9
        solver.parameters["abs_tolerance"] = 1e-12
        solver.parameters["max_iter"]      = 25
        solver.parameters["inner_rel_tolerance"] = self.inner_tol
        solver.parameters["c_armijo"] = 1e-4
        solver.parameters["GN_iter"] = 5
        
        self.x_map = solver.solve(a0)
        
        if solver.converged:
            print "\nConverged in ", solver.it, " iterations."
        else:
            print "\nNot Converged"
        
        print "Termination reason: ", solver.termination_reasons[solver.reason]
        print "Final gradient norm: ", solver.final_grad_norm
        print "Final cost: ", solver.final_cost
        
        return solver.converged
    
    def KLE_GaussianPost(self):
        self.model.setPointForHessianEvaluations(self.x_map)
        Hmisfit = ReducedHessian(self.model, self.inner_tol, gauss_newton_approx=False, misfit_only=True)
        k = 50
        p = 20
        size_param = self.x_map[PARAMETER].array().shape[0]
        print "Double Pass Algorithm. Requested eigenvectors: {0}; Oversampling {1}.".format(k,p)
        Omega = np.random.randn(size_param, k+p)
        self.d_misfit, self.U_misfit = doublePassG(Hmisfit, self.model.prior.R, self.model.prior.Rsolver, Omega, k)
        print "d_misfit \n", self.d_misfit, "\n"
        self.d_misfit[self.d_misfit < 0 ] = 0.
        self.laplace_approx_posterior = GaussianLRPosterior(self.model.prior, self.d_misfit, self.U_misfit)
        self.laplace_approx_posterior.mean = self.x_map[PARAMETER].copy()
        Gamma_post = to_dense( MyOperator(self.laplace_approx_posterior.Hlr, self.model.prior.M) )
        #Gamma_pre = to_dense( MyOperator(self.model.prior.Rsolver, self.model.prior.M) )
        M_dense = to_dense(self.model.prior.M)
        self.d_gaussianPost, self.U_gaussianPost = eigh(Gamma_post,M_dense)
        
        print "d_gaussianPost \n", self.d_gaussianPost, "\n"
        self.d_gaussianPost[self.d_gaussianPost < 0 ] = 0.
        
        return size_param
    
    def Eval(self, eta):
        scaling = eta*np.sqrt(self.d_gaussianPost)
        return np.dot(self.U_gaussianPost, scaling) + self.x_map[PARAMETER].array()
    
    def quit(self):
        if self.is_started:
            self.sock.shutdown(socket.SHUT_RDWR)
        
        self.sock.close()

                