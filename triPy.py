# -*- coding: utf-8 -*-
"""
=== triPy ===
=== Version 1.1 ===
Common calculations for triangular geometric data.  Module includes:

triPy.faceNormals	Returns a NumPy array of face normal components when
					passed NumPy arrays of vertex coordinates and
					connectivity data.

triPy.vertexNormals	Returns a NumPy array of vertex normal components when
					passed NumPy arrays of unique vertex coordinates and
					structured connectivity data.

triPy.inCentres		Returns a NumPy array of face centre coordinates when
					passed NumPy arrays of triangular vertex coordinates
					and connectivity data.

triPy.readASCIISTL	Returns a list object composed of the data from an
					ASCII encoded STL file.

triPy.readBinarySTL	Returns a list object composed of the data from a
					binary encoded STL file.

triPy.readSTLFile	Automatically distinguishes between binary and ASCII
					encoding and calls 'readBinarySTL' or 'readASCIISTL'
					as is appropriate.

triPy.stlToNumpy	Calls 'readSTLFile' and outputs two NumPy '.npy'
					format files containing the vertex and face data of the
					interpreted STL file.

triPy.inferConnect	Infers connections between vertices in unstructured
					triangulated meshes.

triPy.triArea		Calculates surface area of triangular faces defined by three
					sets of XYZ coordinates.

triPy.rotxyz		Rotates a set of triangles represented by XYZ coordinates
					about the global X, Y and Z axes (in that order).

triPy.rotx			Rotates a set of triangles represented by XYZ coordinates
					about the global X axis.

triPy.roty			Rotates a set of triangles represented by XYZ coordinates
					about the global Y axis.

triPy.rotz			Rotates a set of triangles represented by XYZ coordinates
					about the global Z axis.

Created on Mon Jul 22 18:33:08 2013
@author: Nathan Donaldson
"""

def faceNormals(v, f, mode='unit'):
	"""
	=== triPy.faceNormals ===
	Calculates face normal vectors for triangulated geometry data.  Returns an
	m-by-3 NumPy array consisting of UVW components for all faces in input arrays.

	=== Inputs ===
	'v'		An m-by-3 NumPy array of XYZ data representing the coordinates
			of the triangle vertices.
	'f'		An (m/3)-by-3 NumPy array of vertex connectivity data for each
			face.
	'mode'	The magnitude of the returned normals ('unit' returns the unit
			vector normal, 'mag' returns the actual magnitude of the
	normal vector as calculated).  The default is 'unit'.

	=== Usage ===
	import triPy
	n_faces = triPy.faceNormals(v, f, mode='unit')

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

	n_faces = normals

	return n_faces

def vertexNormals(v, f, n):
	"""
	=== triPy.vertexNormals ===
	Calculates vertex normal vectors for triangulated geometry data.  Returns an
	m-by-3 NumPy array consisting of UVW components for all vertices in input arrays.

	=== Inputs ===
	'v'			An m-by-3 NumPy array of XYZ data representing the coordinates
				of unique triangle vertices.
	'f'			A NumPy array of vertex connectivity data for each face (not
				original unconnected STL data).
	'n_faces'	A NumPy array of the mesh face normals

	=== Usage ===
	import triPy
	n_verts = triPy.vertexNormals(v, f, n)

	@author: Nathan Donaldson
	"""

	import numpy as np

	v_neighbours_f = list(-1*(np.ones([len(v), 1], dtype=np.int)))
	n_verts = -1*(np.ones([len(v), 3]))

	# List faces which border a given vertex
	# i.e. v_neighbours_f[vertex number] = faces which border vertex
	for count in range(0, len(v)):
		v_neighbours_f[count] = np.append(v_neighbours_f[count], np.where(f==count)[0])
		# Delete placeholder value at start of each array
		v_neighbours_f[count] = np.delete(v_neighbours_f[count], [0])

	for count2 in range(0, len(v)):
		n_verts[count2, 0] = np.mean(n[v_neighbours_f[count2], 0])
		n_verts[count2, 1] = np.mean(n[v_neighbours_f[count2], 1])
		n_verts[count2, 2] = np.mean(n[v_neighbours_f[count2], 2])

	return n_verts

def inCentres(v, f):
	"""
	=== triPy.inCentres ===
	Calculates face centroids for triangulated geometry data.  Returns an
	m-by-3 NumPy array dubbed 'centroids' consisting of XYZ coordinates for all
	faces in input arrays.

	=== Inputs ===
	'v'		An m-by-3 NumPy array of XYZ data representing the coordinates
			of the triangle vertices.
	'f'		An (m/3)-by-3 NumPy array of vertex connectivity data for each
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
	numFaces = 0
	solidName = list()
	#vertices = np.array([])
	triangles = np.array([])

	# Access STL file
	inputFile = open(filepath)
	print 'Importing ASCII STL file'

	# Read data in file and store as string
	inputFile.seek(0)
	inputStr = inputFile.read()

	# Count number of faces in file for array pre-allocation
	# Iterate through solids (strip 'solid' string)
	for solidStr in re.findall(r'solid\s(.*?)endsolid', inputStr, re.S):
		solidName.append(re.match(r'^(.*)$', solidStr, re.M).group(0))
		print 'Checking object %s' % solidName[object_count]

		# Iterate through facets (strip 'facet' string, normal is unused)
		for facetStr in re.findall(r'facet\s(.*?)endfacet', solidStr, re.S):
			numFaces += 1

	vertices = np.zeros([3 * numFaces, 3])

	# Convert data from string to numpy arrays (strip human readable sections first)
	# Iterate through solids (strip 'solid' string)
	for solidStr in re.findall(r'solid\s(.*?)endsolid', inputStr, re.S):
		solidName.append(re.match(r'^(.*)$', solidStr, re.M).group(0))
		print 'Importing object %s with %i faces' % (solidName[object_count], numFaces)

		# Iterate through facets (strip 'facet' string, normal is unused)
		for face_index, facetStr in enumerate(re.findall(r'facet\s(.*?)endfacet', solidStr, re.S)):

			# Iterate through outer loops (strip 'outer loop' string)
			for outerLoopStr in re.findall(r'outer\sloop(.*?)endloop', facetStr, re.S):

				# Iterate through vertices (strip 'vertex' string, leaving raw vertex data)
				for index, vertexStr in enumerate(re.findall(r'vertex\s(.*)$', outerLoopStr, re.M)):
					#np.append(vertices, [float(coord) for coord in vertexStr.split()])
					vertices[3 * face_index + index, 0] = float(vertexStr.split()[0])
					vertices[3 * face_index + index, 1] = float(vertexStr.split()[1])
					vertices[3 * face_index + index, 2] = float(vertexStr.split()[2])

	#numFaces = len(vertices)/9
	elapsed = time.time() - t
	print 'Import complete - %i faces processed in %f seconds' % (numFaces, elapsed)

	# Reorganise vertex data into m-by-3 array
	#vertices = np.reshape(vertices,(len(vertices)/3,3))

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

	#vertices = np.array([])
	triangles = np.array([])

	# Access STL file
	inputFile = open(filepath, 'rb')
	print 'Importing binary STL file'

	# Read file header (80 bytes, ASCII)
	inputFile.seek(0, 0)
	header = inputFile.read(80)

	# Read number of faces (1 x UINT32, 4 bytes)
	numFaces = struct.unpack('<I', inputFile.read(4))[0]

	vertices = np.zeros([numFaces])

	# Read vertex data
	# Iterate through faces
	for n in range(numFaces):

		# Skip normal vector (3 x REAL32, 12 bytes)
		inputFile.seek(12, 1)

		# Convert vertices to NumPy arrays (9 x REAL32, 36 bytes)
		#vertices = np.append(vertices, struct.unpack('<9f', inputFile.read(36)))
		vertices[n] = struct.unpack('<9f', inputFile.read(36))

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

	'filename_vertices.npy':	A NumPy array file of the XYZ coordinates of the
								geometry's vertices (FLOAT)
	'filename_faces.npy':		A NumPy array file of vertex connectivity data
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

def inferConnect(v, f, n_faces, c):
	"""
	=== triPy.inferConnect ===
	Infers connectivity between points in unstructured triangular grids and
	also calculates vertex normals. Returns the list 'meshConnections' which
	consists of the following 4 cells:

	[0]: A NumPy array of the XYZ coordinates of unique vertices belonging to
		the input mesh.
	[1]: A NumPy array of unique mesh edges, with no duplication between faces.
		The array is structured such that the XYZ coordinates of an edge's
		first point are listed, followed by the XYZ coordinate of the second.
	[2]: A NumPy array of the face connectivity data that relates to the
		vertex data in the first cell i.e. for each triangle, the indexes of
		whichever three vertices are used to construct it.
	[3]: A NumPy array of the mesh vertex normals.

	=== Inputs ===
	'v'			A NumPy array of XYZ data representing the coordinates of the
				unordered triangle vertices.
	'f'			A NumPy array of vertex connectivity data for each face.
	'n_faces'	A NumPy array of the components of the face normals.
	'c'			A NumPy array of the XYZ coordinates of the face centres.

	=== Usage ===
	import triPy
	meshConnection = triPy.inferConnect(v, f, n_faces, c)
	"""

	import numpy as np

	# Find input array lengths
	numVert = len(v)
	numFace = len(f)

	# Construct data structure containing vertex triplets which may be
	# interrogated for unique vertices
	ncols = v.shape[1]
	dtype = v.dtype.descr * ncols
	struct = v.view(dtype)

	# Require a list of vertices with a note of which triangles they belong to
	v_uniq = np.unique(struct)
	tris = list(np.zeros([len(v_uniq), 1]))

	# Find indexes of elements in struct where entries of v_uniq occur
	# Returns list of arrays of face indices which utilise numbered veritces
	# (vertex number is denoted by index of tris list)
	for i in range(0, len(v_uniq)):
		tris[i] = np.where(v_uniq[i]==struct)[0]

	f_new = -1*(np.ones([numVert], dtype=np.int))

	# Re-broadcast face-vertex indexing so that an array of unique vertices referenced
	# by a given face is returned (face number is denoted by index of f_new array)
	for j in range(0, len(tris)):
		for k in tris[j]:
			f_new[k] = j

	f_new = np.reshape(f_new, [numFace, 3])
	v_uniq = v_uniq.view(v.dtype).reshape(-1, ncols)
	v_uniq = np.float64(v_uniq)

	v_neighbours_f = list(-1*(np.ones([len(tris), 1], dtype=np.int)))
	n_verts = -1*(np.ones([len(tris), 3]))

	# List faces which border a given vertex
	# i.e. v_neighbours_f[vertex number] = faces which border vertex
	for l in range(0, len(tris)):
		v_neighbours_f[l] = np.append(v_neighbours_f[l], np.where(f_new==l)[0])
		# Delete placeholder value at start of each array
		v_neighbours_f[l] = np.delete(v_neighbours_f[l], [0])

	#for l in range(0, len(tris)):
	#    v_neighbours_f[l] = np.delete(v_neighbours_f[l], [0])

	v_neighbours_v = list(-1*(np.ones([len(tris), 1], dtype=np.int)))

	# Find vertices which belong to faces bordering given vertex
	# (calculates 'one-ring' neighbourhood around given vertex)
	for m in range(0, len(tris)):
		v_neighbours_v[m] = np.append(v_neighbours_v[m], f_new[v_neighbours_f[m]])
		v_neighbours_v[m] = np.delete(v_neighbours_v[m], [0])
		v_neighbours_v[m] = np.unique(v_neighbours_v[m])
		# Delete centre vertex from list of bordering vertices
		v_neighbours_v[m] = np.delete(v_neighbours_v[m], np.where(v_neighbours_v[m]==m)[0])

	# Calculate vertex normals (mean of surrounding face normals)
	for count2 in range(0, len(tris)):
		n_verts[count2, 0] = np.mean(n_faces[v_neighbours_f[count2], 0])
		n_verts[count2, 1] = np.mean(n_faces[v_neighbours_f[count2], 1])
		n_verts[count2, 2] = np.mean(n_faces[v_neighbours_f[count2], 2])

	index_size = np.zeros([len(v_neighbours_v)], dtype=np.int)

	# Find size of each vertex's one-ring neighbourhood
	# (number of surrounding vertices)
	for count5 in range(0, len(v_neighbours_v)):
		index_size[count5] = len(v_neighbours_v[count5])
	index_size_total = np.int(np.cumsum(index_size)[-1])

	edges = np.zeros([index_size_total, 6])
	index_counter = 0

	index_size_cumsum = np.cumsum(index_size)-index_size[0]

	# Calculate individual mesh edges
	# Edges are categorised by XYZ coordinates of first point, followed by XYZ
	# coordinates of second point
	for count3 in range(0, len(v_neighbours_v)):
		for count4 in range(0, len(v_neighbours_v[count3])):
			edges[index_counter, 0:3] = v_uniq[count3, :]
			edges[index_counter, 3:7] = v_uniq[v_neighbours_v[count3][count4]]
			index_counter = index_counter + 1

	meshConnections = list([v_uniq, f_new, edges, n_verts])

	return meshConnections

def triArea(v, f):
	"""
	=== triPy.triArea ===
	Calculates face area for triangular three dimensional elements.  Returns an
	m-by-1 NumPy array consisting of element areas.  Makes use of NumPy cross
	product function "numpy.cross".

	=== Inputs ===
	'v'		An m-by-3 NumPy array of XYZ data representing the coordinates
			of the triangle vertices.
	'f'		An (m/3)-by-3 NumPy array of vertex connectivity data for each
			face.

	=== Usage ===
	import triPy
	A_faces = triPy.triArea(v, f)

	@author: Nathan Donaldson
	"""

	import numpy as np

	if np.shape(f) == (3, ):
		vector_1 = v[f][2] - v[f][0]
		vector_2 = v[f][1] - v[f][0]
		A = 0.5 * np.linalg.norm(np.cross(vector_1, vector_2))

	else:
		A = np.zeros(len(f))
		vector_1 = np.zeros([len(f), 3])
		vector_2 = np.zeros([len(f), 3])
		for index, value in enumerate(f):
			vector_1[index, :] = v[f[index, :]][2] - v[f[index, :]][0]
			vector_2[index, :] = v[f[index, :]][1] - v[f[index, :]][0]
			A[index] = 0.5 * np.linalg.norm(np.cross(vector_1[index], vector_2[index]))

	return A
	
def rotxyz(points, beta, theta, alpha, mode='deg'):
	"""
	=== triPy.rotxyz ===
	Rotates points about the global X, Y and Z axes by the angles beta, theta, 
	and alpha (roll, yaw, and pitch, respectively) in that order.  Rotation 
	matrices are sourced from http://mathworld.wolfram.com/RotationMatrix.html
	
	=== Inputs ===
	'points'	XYZ coordinates of the points to be rotated
	'beta'		Angle of rotation about the X axis (roll)
	'theta'		Angle of rotation about the Y axis (yaw)
	'alpha'		Angle of rotation about the Z axis (pitch)
	'mode'		Defines whether angles are expressed in radians or degrees 
				(degrees are the default)
	
	=== Usage ===
	import triPy
	rot_x, rot_y, rot_z = rotMatrix(x, y, z, beta, theta, alpha)
	"""
	
	import numpy as np
	
	if mode == 'rad':
		pass
	elif mode == 'deg':
		beta	=	np.deg2rad(beta)
		theta	=	np.deg2rad(theta)
		alpha	=	np.deg2rad(alpha)
	else:
		print 'ERROR: Incorrect angle type specified.  Assuming degrees.'
	
	# Rotation about X axis
	x_rot_mat = np.array([	[1, 	0, 				0			],
							[0, 	np.cos(beta), 	np.sin(beta)],
							[0, 	-np.sin(beta), 	np.cos(beta)]	])
	
	# Rotation about Y axis
	y_rot_mat = np.array([	[np.cos(theta), 0, 	-np.sin(theta)	],
							[0, 			1,	0				],
							[np.sin(theta), 0, 	np.cos(theta)	]	])
	
	# Rotation about Z axis
	z_rot_mat = np.array([	[np.cos(alpha),		np.sin(alpha),	0],
							[-np.sin(alpha),	np.cos(alpha),	0],
							[0,					0,				1]	])
	
	# Sequentially rotate input points about X, Y and then Z axes
	if np.size(points) == 3:
		rot_x 	= np.dot(points, x_rot_mat)
		rot_xy 	= np.dot(rot_x, y_rot_mat)
		rot_xyz = np.dot(rot_xy, z_rot_mat)
	else:
		rot_x 	= np.dot(points, x_rot_mat)
		rot_xy 	= np.dot(rot_x, y_rot_mat)
		rot_xyz = np.dot(rot_xy, z_rot_mat)
		
	return rot_xyz
	
def rotx(points, beta, mode='deg'):
	"""
	=== triPy.rotx ===
	Rotates points about the global X axis by the angles beta (roll).  Rotation 
	matrices are sourced from http://mathworld.wolfram.com/RotationMatrix.html
	
	=== Inputs ===
	'points'	XYZ coordinates of the points to be rotated
	'beta'		Angle of rotation about the X axis (roll)
	'mode'		Defines whether angles are expressed in radians or degrees 
				(degrees are the default)
	
	=== Usage ===
	import triPy
	rot_x = rotx(points, beta)
	"""
	
	import numpy as np
	
	if mode == 'rad':
		pass
	elif mode == 'deg':
		beta	=	np.deg2rad(beta)
	else:
		print 'ERROR: Incorrect angle type specified.  Assuming degrees.'
	
	# Rotation about X axis
	x_rot_mat = np.array([	[1, 	0, 				0			],
							[0, 	np.cos(beta), 	np.sin(beta)],
							[0, 	-np.sin(beta), 	np.cos(beta)]	])
	
	# Sequentially rotate input points about X, Y and then Z axes
	if np.size(points) == 3:
		rot_x 	= np.dot(points, x_rot_mat)
	else:
		rot_x 	= np.dot(points, x_rot_mat)
		
	return rot_x
	
def roty(points, theta, mode='deg'):
	"""
	=== triPy.roty ===
	Rotates points about the global Y axis by the angles theta (yaw).  Rotation 
	matrices are sourced from http://mathworld.wolfram.com/RotationMatrix.html
	
	=== Inputs ===
	'points'	XYZ coordinates of the points to be rotated
	'theta'		Angle of rotation about the Y axis (yaw)
	'mode'		Defines whether angles are expressed in radians or degrees 
				(degrees are the default)
	
	=== Usage ===
	import triPy
	rot_y = roty(points, theta)
	"""
	
	import numpy as np
	
	if mode == 'rad':
		pass
	elif mode == 'deg':
		theta	=	np.deg2rad(theta)
	else:
		print 'ERROR: Incorrect angle type specified.  Assuming degrees.'
	
	# Rotation about Y axis
	y_rot_mat = np.array([	[np.cos(theta), 0, 	-np.sin(theta)	],
							[0, 			1,	0				],
							[np.sin(theta), 0, 	np.cos(theta)	]	])
	
	# Sequentially rotate input points about X, Y and then Z axes
	if np.size(points) == 3:
		rot_y 	= np.dot(points, y_rot_mat)
	else:
		rot_y 	= np.dot(points, y_rot_mat)
		
	return rot_y
	
def rotz(points, alpha, mode='deg'):
	"""
	=== triPy.rotz ===
	Rotates points about the global Z axis by the angles theta (pitch).  Rotation 
	matrices are sourced from http://mathworld.wolfram.com/RotationMatrix.html
	
	=== Inputs ===
	'points'	XYZ coordinates of the points to be rotated
	'alpha'		Angle of rotation about the Z axis (pitch)
	'mode'		Defines whether angles are expressed in radians or degrees 
				(degrees are the default)
	
	=== Usage ===
	import triPy
	rot_z = rotz(points, alpha)
	"""
	
	import numpy as np
	
	if mode == 'rad':
		pass
	elif mode == 'deg':
		alpha	=	np.deg2rad(alpha)
	else:
		print 'ERROR: Incorrect angle type specified.  Assuming degrees.'
	
	# Rotation about Z axis
	z_rot_mat = np.array([	[np.cos(alpha),		np.sin(alpha),	0],
							[-np.sin(alpha),	np.cos(alpha),	0],
							[0,					0,				1]	])
	
	# Sequentially rotate input points about X, Y and then Z axes
	if np.size(points) == 3:
		rot_z 	= np.dot(points, z_rot_mat)
	else:
		rot_z 	= np.dot(points, z_rot_mat)
		
	return rot_z
