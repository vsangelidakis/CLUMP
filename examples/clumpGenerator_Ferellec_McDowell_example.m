%% Script to run clumpGenerator_Ferellec_McDowell
% Copyright © 2020 V. Angelidakis.
clc; clear; close all
addpath(genpath('../functions'))  % Add path to dependencies (external codes)

%% Load particle shape from stl
% stlFile='Octahedron_Coarse_Mesh.stl'; dmin=0.01;	rmin=0.01;	rstep=0.001;	pmax=1.0;	seed=5;	output='FM_octaCoarse.txt';	visualise=true;
stlFile='Octahedron_Fine_Mesh.stl';     dmin=0.01;	rmin=0.01;	rstep=0.001;	pmax=1.0;	seed=5;	output='FM_octaFine.txt';	visualise=true;
% stlFile='Hexahedron_Coarse_Mesh.stl'; dmin=0.01;	rmin=0.01;	rstep=0.001;	pmax=1.0;	seed=5;	output='FM_hexaCoarse.txt';	visualise=true;
% stlFile='Hexahedron_Fine_Mesh.stl';   dmin=0.01;	rmin=0.01;	rstep=0.001;	pmax=1.0;	seed=5;	output='FM_hexaFine.txt';	visualise=true;

tic
[mesh, clump]=clumpGenerator_Ferellec_McDowell( stlFile, dmin, rmin, rstep, pmax, seed, output, visualise );
toc

disp(['Total number of spheres: ', num2str(clump.numSpheres)])