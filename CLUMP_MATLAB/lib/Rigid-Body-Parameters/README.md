# Rigid Body Parameters

[![View Rigid body parameters of closed surface meshes on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/48913-rigid-body-parameters-of-closed-surface-meshes)

In order to simulate dynamic behaviour of a rigid-body, one requires knowledge of a set of rigid-body parameters
such as the total mass of the rigid-body, the center of mass, as well as the moments and products of inertia.
The purpose of this submission is to provide a function which computes exact rigid-body parameters of objects 
represented by closed, triangular surface meshes. The principles underlying the calculations are based on the 
divergence theorem and are explained in detail in the attached .pdf document. This submission also includes two 
functions that take as input an arbitrary mesh and output parameters of a primitive object, such as an ellipsoid
or a cuboid, with exactly the same inertial parameters as the input object. Finally, `VisualizeLocaFrame.m` 
function can be used for visualizing local frames of reference constructed from principal axes of inertia.

## Quick Demo

	load('sample_mesh') 	
	RBP=RigidBodyParams(TR);
	disp(RBP)
	VisualizeLocalFrame(TR)

## Enforcing Consistent and Proper Face Orientation

All calculations are based on the assumption that the input mesh is closed, manifold, and has outward pointing normals.
To obtain outward point normals, the vertices of all faces must have counterclockwise ordering. If you know or suspect
the input the mesh has either inconsistent or improper face orientation use function `ConsistentNormalOrientation` prior
to computing the rigid-body parameters. Here is an example:


	load('sample_mesh')
	[F,V]=GetMeshData(TR);
	
	% Randomly mix-up orientations of the faces to simulate the problem above
	Nf=size(F,1);
	idx=randn(Nf,1)>0;	
	F2=F;
	F2(idx,:)=fliplr(F(idx,:));
	TR2=triangulation(F2,V);
	fprintf('\nNumber of inverted faces: %d\n',nnz(idx))
		
	% Enforce proper face orientation
	[TR2_fix,cnt]=ConsistentNormalOrientation(TR2);
	fprintf('Number of faces corrected: %d\n\n',cnt)

	% Verify that the output is identical to RBP for the original mesh
	RBP_fix=RigidBodyParams(TR2_fix);
	
	fprintf('Corrected mesh:\n')
	disp(RBP_fix)
	
	fprintf('Original (reference) mesh:\n')
	disp(RigidBodyParams(TR))

	% Result you would have gotten without ensuring proper face orientation:
	RBP2=RigidBodyParams(TR2);
	fprintf('Uncorrected mesh:\n')
	disp(RBP2)
 


## License
[MIT] © 2019 Anton Semechko 
a.semechko@gmail.com

[MIT]: https://github.com/AntonSemechko/Rigid-Body-Parameters/blob/master/LICENSE.md
