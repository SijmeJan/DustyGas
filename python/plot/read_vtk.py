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

filename = '../../data/state-10-rank-0.vtk'

s = SnapShot(filename)
cs = plt.tricontourf(s.x[:,0], s.x[:,1], s.Q[:,0], 100)
plt.colorbar(cs)

#n = 21
#e = np.zeros((n))

#for i in range(0, n):
#    s = SnapShot('../../data/state-{}-rank-0.vtk'.format(i))
#    e[i] = np.sum(s.Q[:,6]*s.Q[:,6]/s.Q[:,4])

#plt.yscale('log')
#plt.plot(e)

#plt.plot(0.1*np.exp(0.21*np.arange(0,n)))


plt.show()
