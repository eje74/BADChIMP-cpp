import numpy as np
from struct import pack,unpack 
import sys

###
###
###
def read_data(f, dtype, stride, N):
    unpack_code = {'unsigned_char':'B', 'char':'b', 'short':'h', 'int':'i',
                   'long':'q', 'double':'d', 'float':'f'}
    data_length = {'unsigned_char':1, 'char':1, 'short':2, 'int':4,
                   'long':8, 'double':8, 'float':4}
    elements = N*stride
    size = elements*data_length[dtype]
    data = f.read(size)

    ### Check length of data read, and abort if sizes don't match
    if len(data) < size:
        print('WARNING! Unmatched size for ' + filename.split('/')[-1] + ', file skipped!')  
        return(-1)
    
    data = np.asarray(unpack('>'+str(elements)+unpack_code[dtype], data))
    if stride > 1:
        data = data.reshape([stride, int(N)], order='F').transpose(1,0)

    data
    return(data)

###
###
###
def read(filename, rescale=1e-3):
    ### open file
    #print 'Opening %s' % filename    
    f = open(filename, 'rb')

    vtk = {}

    ### read ASCII header
    f.readline() # # vtk DataFile Version
    f.readline() # title
    f.readline() # BINARY
    f.readline() # DATASET POLYDATA

    ### read POINTS
    (type, N, dtype) = f.readline().split() # POINTS num float
    stride = 3
    vtk['points'] = read_data(f, dtype.decode(), stride, int(N))
    if rescale > 0:
        vtk['points'] *= rescale
    
    ### read LINES
    dtype = 'int'
    stride = 1
    f.readline() # newline
    (type, nlines, N) = f.readline().split() # LINES numlines numdata
    lines = read_data(f, dtype, stride, int(N))

    ### save as separate lines
    vtk['line'] = []
    vtk['length'] = []
    pos = 0
    for i in range(int(nlines)):
        n = lines[pos]
        pos += 1
        vtk['line'].insert(i, lines[pos:pos+n])
        vtk['length'].insert(i, np.zeros(n))
        vtk['length'][i][1:] = np.cumsum(np.sqrt(np.sum(np.diff(vtk['points'][ vtk['line'][i] ], axis=0)**2, axis=1))) 
        #vtk['length'].insert(i, np.cumsum(np.sqrt(np.sum(np.diff(vtk['points'][ vtk['line'][i],: ], axis=0)**2, axis=1))) )
        pos += n
        
    vtk['line'] = np.asarray(vtk['line'])
    ### calculate length of 
        
    ### read FIELD
    stride = 1
    f.readline()  # newline
    f.readline()  # POINT_DATA num
    f.readline()  # FIELD FieldData 1
    (name, nfields, N, dtype) = f.readline().split() # name num elements dtype
    vtk[name.decode()] = read_data(f, dtype.decode(), stride, int(N))
    
    f.close()
    return(vtk)


###
###
###
def append_fielddata(filename, array, name):
    with open(filename, 'rb') as f:
        s = f.read()
    ### find position of 'FieldData' heading
    pattern = 'FieldData '
    p = s.find(str.encode(pattern))
    if p<0:
        print('"%s" not found in %s, aborting...' % (pattern, filename))
        f.close()
        return(-1)
    ### read number
    n = int(s[p+10:p+11])
    ### increment and replace number
    s1 = '%s%d' % (pattern,n)
    s2 = '%s%d' % (pattern,n+1)
    ss = s.replace(str.encode(s1), str.encode(s2))
    ss += str.encode('%s 1 %d double\n' % (name, len(array)))
    ss += pack('>%s%s' % (len(array), 'd'), *np.reshape(array, len(array), order='F'))
    newfile = filename.split('.')[0]+'_%s.vtk' % name
    with open(newfile, 'wb') as f:
        f.write(ss)
    print('   Wrote %s' % newfile)
    return(1)
    #f.write(str.encode('Pressure 1 %d double\n' % len(array)))
    #f.write(pack('>%s%s' % (len(array), 'd'), *reshape(array, len(array), order='F')))
    #f.close()
    
    
    

#if __name__ == "__main__":
#    read_centerlines(sys.argv[1])
