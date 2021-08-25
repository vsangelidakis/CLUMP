%% Script to run clumpGenerator_Euclidean_3D
% 2021 © V. Angelidakis, S. Nadimi, M. Otsubo, S. Utili.
clc; clear; close all
addpath(genpath('../functions'))  % Add path to in-house functions

%% Load particle shape from stl file
inputGeom='ParticleGeometries/Octahedron_Coarse_Mesh.stl';   N=24; rMin=0; div=102; overlap=0.6; output='EU_octaCoarse.txt'; visualise=true;
% inputGeom='ParticleGeometries/Hexahedron_Coarse_Mesh.stl'; N=24; rMin=0; div=102; overlap=0.6; output='EU_hexaCoarse.txt'; visualise=true;

%% Load particle shape from binary 3D image
% inputGeom='ParticleGeometries/matFile.mat';   N=24; rMin=0; div=102; overlap=0.6; output='FM_3D_image.txt'; visualise=true;

tic
[mesh, clump]=GenerateClump_Euclidean_3D( inputGeom, N, rMin, div, overlap, output, visualise ); % true to plot mesh and clump
toc

disp(['Total number of spheres: ', num2str(clump.numSpheres)])
