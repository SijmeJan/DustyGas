import numpy as np
import os

from vtk import vtkUnstructuredGridReader, vtkCellCenters
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
from matplotlib import pyplot as plt
from scipy.special import roots_legendre

class SnapShot():
    def __init__(self, filename):
        reader = vtkUnstructuredGridReader()
        reader.SetFileName(filename)
        reader.ReadAllVectorsOn()
        reader.ReadAllScalarsOn()
        reader.Update()

        output = reader.GetOutput()

        # State vector
        self.Q = vtk_to_numpy(output.GetCellData().GetArray(0))

        # Cell centre coordinates
        x = np.zeros((np.shape(self.Q)[0],3)) - 1
        temp = numpy_to_vtk(x)
        vtkCellCenters().ComputeCellCenters(output, temp)
        self.x = vtk_to_numpy(temp)

    def remove_ghost(self, order):
        x_sort = np.sort(np.unique(self.x[:,0]))
        y_sort = np.sort(np.unique(self.x[:,1]))

        sel = np.asarray((self.x[:,0] > x_sort[order-1])*(self.x[:,0] < x_sort[-order])*(self.x[:,1] > y_sort[order-1])*(self.x[:,1] < y_sort[-order])).nonzero()

        self.Q = self.Q[sel]
        self.x = self.x[sel]

class Project():
    def __init__(self, direc, order=2, n_ghost=2):
        self.direc = direc
        self.vtk_files = \
          [x for x in os.listdir(direc) if x.endswith('.vtk')]
        self.base_name = self.vtk_files[0].split('-')[0]
        self.order = order
        self.n_ghost = n_ghost

    def file_name(self, n):
        return self.base_name + '-' + str(n) + '.vtk'

    def file_exists(self, n):
        return self.file_name(n) in self.vtk_files

    def vars(self, s):
        n_vars = len(s.Q[0,:])
        ret = []
        for i in range(0, n_vars):
            ret.append(s.Q[:, i])

        return ret

    def time(self, dt):
        n_file = len(self.vtk_files)
        return np.linspace(0, dt*n_file, n_file)

    def contour(self, n, ax, comp=4):
        if self.file_exists(n):
            s = SnapShot(self.direc + '/' + self.file_name(n))
            s.remove_ghost(self.order*self.n_ghost)

            f = None
            if comp % 4 == 0:
                f = s.Q[:,comp]
            else:
                f = s.Q[:,comp]/s.Q[:,comp - comp % 4]
            return ax.tricontourf(s.x[:,0], s.x[:,1], f, 100)
        else:
            print('File {} not found!'.format(self.file_name(n)))

    def time_evol(self, func, dt, ax):
        n_file = len(self.vtk_files)
        t = self.time(dt)
        e = 0.0*t

        for n in range(0, n_file):
            s = SnapShot(self.direc + '/' + self.file_name(n))
            s.remove_ghost(self.order*self.n_ghost)

            e[n] = func(self.vars(s))

        return ax.plot(t, e)

#order = 2
#n_ghost = 2

#p = Project('../../data/order2_limited_muscl_DMP-4-3',
#            order=2, n_ghost=2)

#fig = plt.figure()
#ax = fig.add_subplot(1, 1, 1)

#ax.set_yscale('log')

# Energy associated with z motions
#f = lambda x: np.sum(x[6]*x[6]/x[4])/len(x[4])
#p.time_evol(f, 0.1, ax)

#cs = p.contour(260, ax)
#plt.colorbar(cs)

#plt.show()

#exit(0)

#sel = np.asarray(s.x[:,1] == np.min(s.x[:,1])).nonzero()
#plt.plot(s.x[sel,0], s.Q[sel,4], marker='o', linestyle='None', color='blue')
#plt.plot(s.x[sel,0], s.Q[sel,3], marker='o', linestyle='None', color='blue')
#plt.plot(s.x[sel,0], s.Q[sel,2], marker='o', linestyle='None', color='green')
#plt.plot(s.x[sel,0] + 0.010471975511966, s.Q[sel,4], marker='x', linestyle='None', color='green')

#n_dust = len(s.Q[0,:])/4 - 1

#if n_dust > 1:
#    Stokes_range = [0.0001, 0.1]
    # Get nodes and weights for Gauss-Legendre
#    x, w = roots_legendre(n_dust)
    # Adjust for integrals from 0 to 1
#    lk = 0.5*(x + 1)

    # Work with log(St)
#    s_range = np.log(np.asarray(Stokes_range))

#    Stokes = np.exp(s_range[0] + (s_range[1] - s_range[0])*lk)

#    plt.xscale('log')
#    plt.yscale('log')

#    plt.plot(Stokes, s.Q[0, 4:-1:4])

#plt.show()
#exit(0)

#n = 108
#e = np.zeros((n))
#direcs = ['../../data/order2', '../../data/order3']
#direcs = ['../../data/order2_limited_muscl_DMP-6-10',
#          '../../data/order2_limited_godunov_DMP-6-10']

#t = np.linspace(0, 0.1*n, n)

#plt.xlabel(r'$\Omega t$')
#plt.ylabel(r'$\int \rho_{\rm d} v_{z,{\rm d}}^2$')
#plt.title('Monodisperse linA')

#for direc in direcs:
#    for i in range(0, n):
#        s = SnapShot(direc + '/state-{}.vtk'.format(i))
#        s.remove_ghost(order*n_ghost)

#        e[i] = np.sum(s.Q[:,6]*s.Q[:,6]/s.Q[:,4])/len(s.Q[:,4])
        #e[i] = np.mean(s.Q[:,1]/s.Q[:,0])

#    plt.yscale('log')
#    plt.plot(t, e)

#plt.plot(t, 0.0000000001*np.exp(0.42*2*t))

#plt.show()
