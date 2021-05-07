import numpy as np
import os

from vtk import vtkUnstructuredGridReader, vtkCellCenters
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
from matplotlib import pyplot as plt
from scipy.special import roots_legendre

class SnapShot():
    def __init__(self, filename):
        print('Reading {}'.format(filename))

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
        #print('Number of cells: ', len(x_sort), len(y_sort))

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
        return dt*np.arange(0, n_file)

    def contour(self, n, ax, comp=4):
        if self.file_exists(n):
            s = SnapShot(self.direc + '/' + self.file_name(n))
            s.remove_ghost(self.order*self.n_ghost)

            f = None

            if comp % 4 == 0:
                f = s.Q[:,comp]
                print(f)
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

    def slice(self, n, ax, direction=0, comp=4, indx=None):
        if self.file_exists(n):
            s = SnapShot(self.direc + '/' + self.file_name(n))
            s.remove_ghost(self.order*self.n_ghost)

            f = None

            if comp % 4 == 0:
                f = s.Q[:,comp]
            else:
                f = s.Q[:,comp]/s.Q[:,comp - comp % 4]

            print(np.shape(s.x))

            return ax.plot(s.x[:,direction], f, marker='o',linestyle='None')
        else:
            print('File {} not found!'.format(self.file_name(n)))
