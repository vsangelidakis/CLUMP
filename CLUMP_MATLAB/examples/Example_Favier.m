%% Script to run GenerateClump_Favier
% 2021 © V. Angelidakis, S. Nadimi, M. Otsubo, S. Utili.
clc; clear; close all
addpath(genpath('../functions'))  % Add path to dependencies (external codes)
addpath(genpath('../lib'))  % Add path to dependencies (external codes)

%% Option 1: Load particle shape from stl file
inputGeom='ParticleGeometries/Ellipsoid_R_2.0_1.0_1.0.stl';
N=10;
chooseDistance='min';
output='FA_Ellipsoid_2.0_1.0_1.0.txt';
visualise=true;

%% Option 2: Load particle shape from binary 3D image
% inputGeom='ParticleGeometries/Binary_3D_Cuboid.mat';
% N=10; 
% chooseDistance='avg';
% output='FA_3D_image.txt';
% visualise=true;

tic 
[mesh, clump]=GenerateClump_Favier( inputGeom, N, chooseDistance, output, visualise );
toc

disp(['Total number of spheres: ', num2str(clump.numSpheres)])
