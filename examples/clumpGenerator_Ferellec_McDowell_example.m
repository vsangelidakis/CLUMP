%% Script to run clumpGenerator_Ferellec_McDowell
% Copyright © 2020 V. Angelidakis. All rights reserved
clc; clear; close all
addpath(genpath('../functions'))  % Add path to dependencies (external codes)

% Uncomment any of the lines below to generate clumps for different particles
stlFile='Octahedron_Coarse_Mesh.stl';	dmin=0.0;		rmin=0.01;	pmax=1.0;	seed=5;	output='octaCoarse.txt';
% stlFile='Octahedron_Fine_Mesh.stl';	dmin=0.0001;	rmin=0.001;	pmax=1.0;	seed=5;	output='octaFine.txt';

% stlFile='Hexahedron_Coarse_Mesh.stl';	dmin=0.0;	rmin=0.01;	pmax=1.0;	seed=5;	output='hexaCoarse.txt';
% stlFile='Hexahedron_Fine_Mesh.stl';	dmin=0.01;	rmin=0.01;	pmax=1.0;	seed=5;	output='hexaFine.txt';

% stlFile='stl_stl_5000_faces.stl'; dmin=2; rmin=20; pmax=1.0; seed=5; output='stl_1000.txt';
% stlFile='LB2_Original_voxel_scale_13078_faces.stl';  dmin=0.34; rmin=5; pmax=1.0; seed=5; output='LB2_1000.txt';
% stlFile='LB3_Original_voxel_scale_9248_faces.stl';  dmin=0.25; rmin=5; pmax=1.0; seed=5; output='LB3_1000.txt';

[mesh, clump]=clumpGenerator_Ferellec_McDowell( stlFile, dmin, rmin, pmax, seed, output );

disp(['Total number of spheres: ', num2str(clump.numSpheres)])

disp('Done!')