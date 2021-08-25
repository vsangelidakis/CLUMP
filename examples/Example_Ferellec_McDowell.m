%% Script to run clumpGenerator_Ferellec_McDowell
% 2021 © V. Angelidakis, S. Nadimi, M. Otsubo, S. Utili.
clc; clear; close all
addpath(genpath('../functions'))  % Add path to in-house functions
addpath(genpath('../lib'))  % Add path to in-house functions

%% Option 1: Load particle shape from stl file
% inputGeom='ParticleGeometries/Octahedron_Coarse_Mesh.stl';	dmin=0.01;	rmin=0.01;	rstep=0.001;	pmax=1.0;	seed=5;	output='FM_octaCoarse.txt';	visualise=true;
% inputGeom='ParticleGeometries/Octahedron_Fine_Mesh.stl';		dmin=0.01;	rmin=0.01;	rstep=0.001;	pmax=1.0;	seed=5;	output='FM_octaFine.txt';	visualise=true;
% inputGeom='ParticleGeometries/Hexahedron_Coarse_Mesh.stl';	dmin=0.01;	rmin=0.01;	rstep=0.001;	pmax=1.0;	seed=5;	output='FM_hexaCoarse.txt';	visualise=true;
% inputGeom='ParticleGeometries/Hexahedron_Fine_Mesh.stl';		dmin=0.01;	rmin=0.01;	rstep=0.001;	pmax=1.0;	seed=5;	output='FM_hexaFine.txt';	visualise=true;
inputGeom='ParticleGeometries/Torus.stl';						dmin=0.1;	rmin=0.01;	rstep=0.01;	pmax=1.0;	seed=5;	output='FM_torus.txt';	visualise=true;

%% Option 2: Load particle shape from binary 3D image
% inputGeom='ParticleGeometries/matFile.mat';	dmin=0.2;	rmin=0.5;	rstep=0.2;	pmax=1.0;	seed=5;	output='FM_hexaFineMoved.txt';	visualise=true;

tic 
[mesh, clump]=GenerateClump_Ferellec_McDowell( inputGeom, dmin, rmin, rstep, pmax, seed, output, visualise );
toc

disp(['Total number of spheres: ', num2str(clump.numSpheres)])