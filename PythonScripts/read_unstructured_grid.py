#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

from vtk import vtkXMLUnstructuredGridReader, vtkPointDataToCellData
from numpy import around, zeros, array
from vtk.util.numpy_support import vtk_to_numpy


def has_array(_name, _data):
    for n in range(_data.GetNumberOfArrays()):
        aname = _data.GetArrayName(n)
        #print(aname, _name)
        if _name in aname:
            return True
    return False



class Grid:
    def __init__(self, _sys_bnd, _ncells, _cell_bnd):
        self.cells = zeros((_ncells, 3))
        self.sys_size = _sys_bnd[1::2]-_sys_bnd[::2]
        self.voxel_size = _cell_bnd[1::2]-_cell_bnd[::2]
        self.dim = around( (_sys_bnd[1::2]-_sys_bnd[::2])/self.voxel_size ).astype('int')
        self.shift = -_sys_bnd[::2]
        self.bounds = zeros(6)
        self.bounds[1::2] = self.sys_size
        print('  Bounds: ', end='')
        print(self.bounds)
        if abs(self.shift).sum() > 0:
            print('  Bounds shifted ', end='')
            print(self.shift, end='')
            print(' from the original bounds ')
            #print(_sys_bnd, end='')
        print('  Dimensions: ', end='')
        print(self.dim)
        
    def get_cell_index(self, _nr):
        return tuple(around(self.cells[_nr]/self.voxel_size).astype('int'))
        
    def set_cell(self, _i, _bound):        
        self.cells[_i] = _bound[::2] + self.shift


def get_geo(file_name):        

    # Read the source file.
    reader = vtkXMLUnstructuredGridReader()
    reader.SetFileName(file_name)
    reader.Update()
    output = reader.GetOutput()

    # name of the array with geometry information
    geoname = 'vtkValidPointMask'

    # check if array is a cell- or point-data array
    cell_data = output.GetCellData()
    cell = has_array(geoname, cell_data)

    # convert to cell-data if array is point-data
    if not cell:
        point_data = output.GetPointData()
        point = has_array(geoname, point_data)
        if not point:
            print('  Array {:s} is missing'.format(geo))
            raise SystemExit('  Aborting...')
        print('  Converting point data to cell data')
        p2c = vtkPointDataToCellData()
        p2c.SetInputData(output)
        #p2c.PassPointDataOn()
        p2c.Update()
        cell_data = p2c.GetOutput().GetCellData()

    # set up grid
    sys_bound = array(output.GetBounds())
    ncells = output.GetNumberOfCells()
    cell_bound = zeros(6)
    output.GetCellBounds(0, cell_bound)
    grid = Grid(sys_bound, ncells, cell_bound)
    for i in range(ncells):
        bound = zeros(6)
        output.GetCellBounds(i, bound)
        grid.set_cell(i, bound)    
        
    # read cell data
    geo = {}
    for n in range(cell_data.GetNumberOfArrays()):
        name = cell_data.GetArrayName(n)
        if geoname not in name:
            continue
        #print(name)
        A = vtk_to_numpy(cell_data.GetArray(n))
        data = zeros(grid.dim, A.dtype)
        for i in range(A.shape[0]):
            data[ grid.get_cell_index(i) ] = A[ i ]
        geo = data

    return geo
        
        
if __name__ == '__main__':
#    import sys
#    from myvtk.IO import write, write_geo_dat
#    from pathlib import Path
    
    # The source file
#    filename = sys.argv[1]
    filename = "/home/ejette/Programs/Python/blateTest3.vtu"

    print()
    geo = get_geo(filename)
    
#    outfile = 'vtkimage_geo.vtk'
#    write({'geo':geo.astype('uint8')}, outfile)
#    print('  Saved {:s}'.format(outfile))
#    print()

#    write_geo_dat(geo, filename=Path(filename).stem+'.dat', resolution=1, echo=True)
