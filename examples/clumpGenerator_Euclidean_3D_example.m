%% Script to run clumpGenerator_Euclidean_3D
% Copyright © 2020 V. Angelidakis.
clc; clear; close all
addpath(genpath('../functions'))  % Add path to dependencies (external codes)

	%% FIXME: The clump is not centered to the same centroid as the particle!! Fix this!!!

%% Load particle shape from stl
stlFile='Octahedron_Coarse_Mesh.stl';   N=24; rMin=0; div=102; overlap=0.6; output='EU_octaCoarse.txt'; visualise=true;
% stlFile='Octahedron_Fine_Mesh.stl';   N=24; rMin=0; div=102; overlap=0.6; output='EU_octaFine.txt'; visualise=true;
% stlFile='Hexahedron_Coarse_Mesh.stl'; N=24; rMin=0; div=102; overlap=0.6; output='EU_hexaCoarse.txt'; visualise=true;
% stlFile='Hexahedron_Fine_Mesh.stl';   N=24; rMin=0; div=102; overlap=0.6; output='EU_hexaFine.txt'; visualise=true;

tic
[mesh, clump]=clumpGenerator_Euclidean_3D( stlFile, N, rMin, div, overlap, output, visualise ); % true to plot mesh and clump
toc

disp(['Total number of spheres: ', num2str(clump.numSpheres)])