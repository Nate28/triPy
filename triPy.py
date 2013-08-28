# -*- coding: utf-8 -*-
"""
=== triPy ===
Common calculations for triangular geometric data.  Module includes:

triPy.normals       -   Returns a NumPy array of face normal components when 
                        passed NumPy arrays of vertex coordinates and 
                        connectivity data
                  
triPy.inCentres     -   Returns a NumPy array of face centre coordinates when 
                        passed NumPy arrays of triangular vertex coordinates 
                        and connectivity data   
                    
triPy.readASCIISTL  -   Returns a list object composed of the data from an 
                        ASCII encoded STL file
                    
triPy.readBinarySTL -   Returns a list object composed of the data from a 
                        binary encoded STL file
                        
triPy.readSTLFile   -   Automatically distinguishes between binary and ASCII
                        encoding and calls 'readBinarySTL' or 'readASCIISTL'
                        as is appropriate
                        
triPy.stlToNumpy    -   Calls 'readSTLFile' and outputs two NumPy '.npy'
                        format files containing the vertex and face data of the 
                        interpreted STL file

Created on Mon Jul 22 18:33:08 2013
@author: Nathan Donaldson
"""

def normals(v, f, mode='unit'):
    """
    === triPy.normals ===    
    Calculates face normal vectors for triangulated geometry data.  Returns an
    m-by-3 NumPy array dubbed 'normals' consisting of UVW components for all 
    faces in input arrays.
    
    === Inputs === 
    'v'       - An m-by-3 NumPy array of XYZ data representing the coordinates 
                of the triangle vertices.
    'f'       - An (m/3)-by-3 NumPy array of vertex connectivity data for each 
                face.
    'mode'    - The magnitude of the returned normals ('unit' returns the unit
                vector normal, 'mag' returns the actual magnitude of the 
                normal vector as calculated).  The default is 'unit'. 
    
    === Usage ===
    import triPy
    normals = triPy.normals(v, f, mode='unit')  
    
    @author: Nathan Donaldson
    """    
    
    import numpy as np
    
    # Generate index arrays of vertices for each face 
    v_0_index = np.reshape(f[:, 0], [len(f)])
    v_1_index = np.reshape(f[:, 1], [len(f)])
    v_2_index = np.reshape(f[:, 2], [len(f)])
    
    # Define two vectors along the sides of each face
    vec_0 = v[v_0_index] - v[v_1_index]
    vec_1 = v[v_1_index] - v[v_2_index]
    
    # Find cross-product of vectors (the face normal)    
    if mode == 'unit':
        normals = np.cross(vec_0, vec_1)    
        mag = np.sqrt((normals[:, 0]**2) + (normals[:, 1]**2) + \
        (normals[:, 2]**2))
        normals[:,0] = normals[:,0]/mag
        normals[:,1] = normals[:,1]/mag
        normals[:,2] = normals[:,2]/mag
        
    elif mode == 'mag':
        normals = np.cross(vec_0, vec_1)
        
    else:
        normals = np.cross(vec_0, vec_1)    
        mag = np.sqrt((normals[:, 0]**2) + (normals[:, 1]**2) + \
        (normals[:, 2]**2))
        normals[:,0] = normals[:,0]/mag
        normals[:,1] = normals[:,1]/mag
        normals[:,2] = normals[:,2]/mag 
        print 'Mode string incorrect; continuing with default \'unit\' mode' 
    return normals

def inCentres(v, f):
    """
    === triPy.inCentres ===    
    Calculates face centroids for triangulated geometry data.  Returns an
    m-by-3 NumPy array dubbed 'centroids' consisting of XYZ coordinates for all 
    faces in input arrays.
    
    === Inputs === 
    'v'       - An m-by-3 NumPy array of XYZ data representing the coordinates 
                of the triangle vertices.
    'f'       - An (m/3)-by-3 NumPy array of vertex connectivity data for each 
                face.
    
    === Usage ===
    import triPy
    centroids = triPy.inCentres(v, f)
    
    @author: Nathan Donaldson
    """    
    
    import numpy as np
    
    # Split vertex coordinates into m-by-3 arrays according to connectivity array 'f'
    x_tri = (v[:, 0])[f] 
    y_tri = (v[:, 1])[f] 
    z_tri = (v[:, 2])[f] 
    
    # Calculate average of each coordinate over each face
    x_avg = (x_tri[:, 0] + x_tri[:, 1] + x_tri[:, 2])/3
    y_avg = (y_tri[:, 0] + y_tri[:, 1] + y_tri[:, 2])/3
    z_avg = (z_tri[:, 0] + z_tri[:, 1] + z_tri[:, 2])/3
    
    # Compile centroid coordinates
    centroids = np.column_stack((x_avg, y_avg, z_avg)) 
    
    return centroids
    
def readASCIISTL(filepath):
    """
    === triPy.readASCIISTL ===
    Reads and converts ASCII STL triangulated mesh data into NumPy arrays.  
    Returns the list 'data' which consists of the following 4 cells:
    
    [0]: A NumPy array of the XYZ coordinates of the geometry's vertices (FLOAT)
    [1]: A NumPy array of vertex connectivity data (INT32)
    [2]: The number of faces in the geometry (INT32)
    [3]: A list of the solid names in the file
    
    === Usage ===
    import triPy
    data = triPy.readASCIISTL(filepath)    
    """
    
    import numpy as np
    import re
    import struct
    import time
    
    t = time.time()
    
    object_count = 0
    solidName = list()
    vertices = np.array([])
    triangles = np.array([])
    
    # Access STL file
    inputFile = open(filepath) 
    print 'Importing ASCII STL file'
    
    # Read data in file and store as string
    inputFile.seek(0)
    inputStr = inputFile.read()
    
    # Convert data from string to numpy arrays (strip human readable sections first)
    # Iterate through solids (strip 'solid' string)
    for solidStr in re.findall(r'solid\s(.*?)endsolid', inputStr, re.S):
        solidName.append(re.match(r'^(.*)$', solidStr, re.M).group(0))
        print 'Processing object %s' % solidName[object_count]
        
        # Iterate through facets (strip 'facet' string, normal is unused)
        for facetStr in re.findall(r'facet\s(.*?)endfacet', solidStr, re.S):
            
            # Iterate through outer loops (strip 'outer loop' string)
            for outerLoopStr in re.findall(r'outer\sloop(.*?)endloop', facetStr, re.S):
                
                # Iterate through vertices (strip 'vertex' string, leaving raw vertex data)
                for vertexStr in re.findall(r'vertex\s(.*)$', outerLoopStr, re.M):                
                    vertices = np.append(vertices,[float(coord) for coord in vertexStr.split()])
    
    numFaces = len(vertices)/9
    elapsed = time.time() - t
    print 'Import complete - %i faces processed in %f seconds' % (numFaces, elapsed)
                    
    # Reorganise vertex data into m-by-3 array
    vertices = np.reshape(vertices,(len(vertices)/3,3))
    
    # Generate face connectivity data
    triangles = np.reshape(np.int32(np.cumsum(np.ones([len(vertices)]))), \
    ((len(vertices)/3),3))-1
    
    # Pack data
    data = []
    data.append((vertices, triangles, numFaces, solidName))    
    
    return data
    
def readBinarySTL(filepath):
    """
    === triPy.readBinarySTL ===
    Reads and converts binary STL triangulated mesh data into NumPy arrays.  
    Returns the list 'data' which consists of the following 4 cells:
    
    [0]: A NumPy array of the XYZ coordinates of the geometry's vertices (FLOAT)
    [1]: A NumPy array of vertex connectivity data (INT32)
    [2]: The number of faces in the geometry (INT32)
    [3]: The file's 80 byte header
    
    === Usage ===
    import triPy
    data = triPy.readBinarySTL(filepath)    
    """    
    
    import numpy as np
    import re
    import struct
    import time
    
    t = time.time()
    
    vertices = np.array([])
    triangles = np.array([])
    
    # Access STL file
    inputFile = open(filepath, 'rb') 
    print 'Importing binary STL file'
    
    # Read file header (80 bytes, ASCII)
    inputFile.seek(0, 0)
    header = inputFile.read(80)
    
    # Read number of faces (1 x UINT32, 4 bytes)
    numFaces = struct.unpack('<I', inputFile.read(4))[0]
    
    # Read vertex data
    # Iterate through faces
    for n in range(numFaces):
        
        # Skip normal vector (3 x REAL32, 12 bytes) 
        inputFile.seek(12, 1)
        
        # Convert vertices to NumPy arrays (9 x REAL32, 36 bytes)    
        vertices = np.append(vertices, struct.unpack('<9f', inputFile.read(36)))
            
        # Skip attribute byte count (1 x UINT16, 2 bytes)
        inputFile.seek(2, 1)
    
    # Reorganise vertex data into numFaces-by-3 array
    vertices = np.reshape(vertices,(len(vertices)/3,3))
            
    # Generate face connectivity data
    triangles = np.reshape(np.int32(np.cumsum(np.ones([len(vertices)]))), \
    ((len(vertices)/3),3))-1
    
    # Print completion message
    elapsed = time.time() - t
    print 'Import complete - %i faces processed in %f seconds' % (numFaces, elapsed)

    # Pack data
    data = []
    data.append((vertices, triangles, numFaces, header))
    
    return data

def readSTLFile(filepath):    
    """
    === triPy.readSTLFile ===
    Reads and converts STL triangulated mesh data into NumPy arrays.  
    Automatically distinguishes between binary and ASCII format STL files and 
    returns the list 'data' which consists of the following 4 cells:
    
    [0]: A NumPy array of the XYZ coordinates of the geometry's vertices (FLOAT)
    [1]: A NumPy array of vertex connectivity data (INT32)
    [2]: The number of faces in the geometry (INT32)
    [3]: File comments (the 80 byte header for binary files or the solid name 
         for ASCII files)
    
    === Usage ===
    import triPy
    data = triPy.readASCIISTL(filepath)    
    """
    # Open file     
    fid = open(filepath)    
    
    # Define associative dictionary for import of different file types
    filetype = {0 : readBinarySTL,
                1 : readASCIISTL}
    
    # Determine file type (assumed binary if first 5 bytes are not human-readable)           
    if fid.read(5) == 'solid':
        data = filetype[1](filepath)
    else:
        data = filetype[0](filepath)
    
    return data
    
def stlToNumpy(filepath, outname='STL'):
    """
    === triPy.stlToNumpy ===
    Reads and converts STL triangulated mesh data into NumPy arrays.  
    Automatically distinguishes between binary and ASCII format STL files and 
    generates two NumPy array files ('.npy' format):
    
    'filename_vertices.npy': A NumPy array file of the XYZ coordinates of the 
                             geometry's vertices (FLOAT)
    'filename_faces.npy':    A NumPy array file of vertex connectivity data 
                             (INT32)
    
    === Usage ===
    import triPy
    data = triPy.stlToNumpy(filepath, outname='NewFileName')

    Where 'filepath' is the path string of the STL file to be interpreted and 
    'outname' is a string representing the common name of the output files.    
    """
    
    import numpy as np    
    
    # Prompt user for output file name
    # outname = input('Enter output file name: ')
    
    # Read STL data and unpack into vertex and face arrays
    data = readSTLFile(filepath)
    v = data[0][0]
    f = data[0][1]
    
    # Define strings for output file names
    outname_v = outname + '_vertices'
    outname_f = outname + '_faces'
    
    # Output data
    np.save(outname_v, v)
    print 'STL vertex data saved to: %s.npy' % (outname_v)
    np.save(outname_f, f)
    print 'STL face data saved to: %s.npy' % (outname_f)
    
    
