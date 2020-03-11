import os

import pandas as pd
import numpy as np
import vtk
from vtk.util import numpy_support as nps

from . import table
from . import version

class VTK: # capitals needed/preferred?
    """
    Manage vtk files.

    tbd
    """

    def __init__(self,geom):
        """tbd."""
        self.geom = geom

    @staticmethod
    def from_rectilinearGrid(grid,size,origin=np.zeros(3)):
        """Check https://blog.kitware.com/ghost-and-blanking-visibility-changes/ for missing data."""
        coordArray = [vtk.vtkDoubleArray(),vtk.vtkDoubleArray(),vtk.vtkDoubleArray()]
        for dim in [0,1,2]:
            for c in np.linspace(origin[dim],origin[dim]+size[dim],grid[dim]):
                coordArray[dim].InsertNextValue(c)

        geom = vtk.vtkRectilinearGrid()
        geom.SetDimensions(*grid)
        geom.SetXCoordinates(coordArray[0])
        geom.SetYCoordinates(coordArray[1])
        geom.SetZCoordinates(coordArray[2])

        return VTK(geom)


    @staticmethod
    def from_unstructuredGrid(nodes,connectivity,cell_type):
        """
        Create an unstructured grid (mesh).

        connectivity: 1 based at the moment
        cell_type: TRIANGLE, 'QUAD', 'TETRA','HEXAHEDRON'

        """
        vtk_nodes = vtk.vtkPoints()
        vtk_nodes.SetData(nps.numpy_to_vtk(nodes))

        cells = vtk.vtkCellArray()
        cells.SetNumberOfCells(connectivity.shape[0])
        T = np.concatenate((np.ones((connectivity.shape[0],1),dtype=np.int64)*connectivity.shape[1],
                           connectivity),axis=1).ravel()
        cells.SetCells(connectivity.shape[0],nps.numpy_to_vtk(T, deep=True, array_type=vtk.VTK_ID_TYPE))

        geom = vtk.vtkUnstructuredGrid()
        geom.SetPoints(vtk_nodes)
        geom.SetCells(eval('vtk.VTK_{}'.format(cell_type.upper())),cells)

        return VTK(geom)


    @staticmethod
    def from_points(nodes,connectivity,cell_type):
        pass


    def write(self,fname):                                              #ToDo: Discuss how to handle consistently filename extensions
        if  (isinstance(self.geom,vtk.vtkRectilinearGrid)):
            writer = vtk.vtkXMLRectilinearGridWriter()
        elif(isinstance(self.geom,vtk.vtkUnstructuredGrid)):
            writer = vtk.vtkXMLUnstructuredGridWriter()
        elif(isinstance(self.geom,vtk.vtkPolyData)):
            writer = vtk.vtkXMLPolyDataWriter()

        writer.SetFileName('{}.{}'.format(os.path.splitext(fname)[0],
                                          writer.GetDefaultFileExtension()))
        writer.SetCompressorTypeToZLib()
        writer.SetDataModeToBinary()
        writer.SetInputData(self.geom)

        writer.Write()


    def add(data,label=None):
        if   isinstance(data,np.ndarray):
            pass
        elif isinstance(data,pd.DataFrame):
            pass
        elif isinstance(data,table):
            pass


    def __repr__(self):
        """ASCII representation of the VTK data."""
        writer = vtk.vtkDataSetWriter()
        writer.SetHeader('DAMASK.VTK v{}'.format(version))
        writer.WriteToOutputStringOn()
        writer.SetInputData(self.geom)
        writer.Write()
        return writer.GetOutputString()