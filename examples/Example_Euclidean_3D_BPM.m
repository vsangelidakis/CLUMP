%% Script to run GenerateClump_Euclidean_3D
% 2021 © V. Angelidakis, S. Nadimi, M. Otsubo, S. Utili.
clc; clear; close all
addpath(genpath('../functions'))  % Add path to in-house functions

%% Load particle shape from stl file
%inputGeom='ParticleGeometries/Octahedron_Coarse_Mesh.stl';	N=200; rMin=0; div=102; overlap=0.0; output='EU_octaCoarse.txt'; visualise=true; rMax_ratio=0.4;
inputGeom='ParticleGeometries/Human_femur.stl';	N=500; rMin=0; div=102; overlap=0.0; output='Human_femur.txt'; visualise=true; rMax_ratio=0.2;

%% Load particle shape from binary 3D image
%inputGeom='ParticleGeometries/Binary_3D_Cube.mat';   N=200; rMin=0; div=102; overlap=0.0; output='image.txt'; visualise=true; rMax_ratio=0.4;

tic
[mesh, clump]=GenerateClump_Euclidean_3D( inputGeom, N, rMin, div, overlap, output, visualise, rMax_ratio); % true to plot mesh and clump
toc

disp(['Total number of spheres: ', num2str(clump.numSpheres)])