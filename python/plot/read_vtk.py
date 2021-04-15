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

filename = '../../data/state-100-rank-0.vtk'
s = SnapShot(filename)
s.remove_ghost(1)

#cs = plt.tricontourf(s.x[:,0], s.x[:,1], s.Q[:,1], 100)
#plt.colorbar(cs)

sel = np.asarray(s.x[:,1] == np.min(s.x[:,1])).nonzero()
plt.plot(s.x[sel,0], s.Q[sel,1], marker='o', linestyle='None', color='blue')
plt.plot(s.x[sel,0], s.Q[sel,2], marker='o', linestyle='None', color='green')


#n = 100
#e = np.zeros((n))
#f = np.zeros((n))

#for i in range(0, n):
#    s = SnapShot('../../data/state-{}-rank-0.vtk'.format(i))
    #e[i] = np.sum(s.Q[:,6]*s.Q[:,6]/s.Q[:,4])
#    e[i] = np.max(s.Q[:,1]/s.Q[:,0])
#    f[i] = np.max(s.Q[:,3]/s.Q[:,0])

#plt.yscale('log')
#plt.plot(np.abs(e))
#plt.plot(np.abs(f))

#plt.plot(0.0001*np.exp(0.42*np.arange(0,n)))


plt.show()
