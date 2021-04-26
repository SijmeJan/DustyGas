import numpy as np
from vtk import vtkUnstructuredGridReader, vtkCellCenters
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
from matplotlib import pyplot as plt

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

order = 2
n_ghost = 2

filename = '../../data/state-220.vtk'
s = SnapShot(filename)
s.remove_ghost(order*n_ghost)

cs = plt.tricontourf(s.x[:,0], s.x[:,1], s.Q[:,4], 100)
plt.colorbar(cs)

plt.show()

exit(0)

#sel = np.asarray(s.x[:,1] == np.min(s.x[:,1])).nonzero()
#plt.plot(s.x[sel,0], s.Q[sel,4], marker='o', linestyle='None', color='blue')
#plt.plot(s.x[sel,0], s.Q[sel,3], marker='o', linestyle='None', color='blue')
#plt.plot(s.x[sel,0], s.Q[sel,2], marker='o', linestyle='None', color='green')
#plt.plot(s.x[sel,0] + 0.010471975511966, s.Q[sel,4], marker='x', linestyle='None', color='green')


#plt.show()
#exit(0)

n = 204
e = np.zeros((n))
#direcs = ['../../data/order2', '../../data/order3']
direcs = ['../../data', '../../data/order2_limited_godunov']

t = np.linspace(0, 0.1*n, n)

plt.xlabel(r'$\Omega t$')
plt.ylabel(r'$\int \rho_{\rm d} v_{z,{\rm d}}^2$')
plt.title('Monodisperse linA')

for direc in direcs:
    for i in range(0, n):
        s = SnapShot(direc + '/state-{}.vtk'.format(i))
        s.remove_ghost(order*n_ghost)

        e[i] = np.sum(s.Q[:,6]*s.Q[:,6]/s.Q[:,4])/len(s.Q[:,4])
        #e[i] = np.mean(s.Q[:,1]/s.Q[:,0])

    plt.yscale('log')
    plt.plot(t, e)

plt.plot(t, 0.0000000001*np.exp(0.42*2*t))

plt.show()
