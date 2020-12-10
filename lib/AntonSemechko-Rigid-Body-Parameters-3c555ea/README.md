# Rigid Body Parameters
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

## License
[MIT] © 2019 Anton Semechko 
a.semechko@gmail.com

[MIT]: https://github.com/AntonSemechko/Rigid-Body-Parameters/blob/master/LICENSE.md