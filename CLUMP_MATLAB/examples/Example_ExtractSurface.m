%% Script to demonstrate the ExtractSurface function
% 2021 © V. Angelidakis, S. Nadimi, M. Otsubo, S. Utili.
clc; clear; close all
addpath(genpath('../functions'))  % Add path to in-house functions
addpath(genpath('../lib'))  % Add path to external functions

%% Input clump

%% Input format 1: Nx4 matrix with columns "x,y,z,r"
% clump=[1,0,0,1.1;
%        2,1,0,1.1]; 

%% Input format 2: struct with fields "positions", "radii"
clump=struct();
clump.positions=[1,0,0;
		 2,1,0;
		 3,0,0;
		 1,0,1];
clump.radii=[1.1;1.1;1.2;1.2];


%% PART 1: CLUMP
%% Tessellate the surface of the clump
N_sphere=200;
N_circle=100;
visualise=true;
[faces,vertices]=ExtractSurface(clump,N_sphere,N_circle,visualise);
view(65,45)

fv=struct();
fv.faces=faces;
fv.vertices=vertices;

%% Write .stl particle with the generated surface mesh (optional)
stlWrite('Clump_Surface.stl',fv) 

%% PART 2: SHAPE
% Characterise the morphology of the surface mesh using SHAPE
% Clone from: https://github.com/vsangelidakis/SHAPE/ to the current
% directory for the code below to work
% 
% addpath(genpath('SHAPE-master/'))			 % Import the dependencies and structure of SHAPE
% load('SHAPE-master/examples/options.mat'); % Add default options of SHAPE
% p=Particle(vertices,faces,[],[],options);	 % Analyse particle morphology
% 
% % Print some geometrical features of the analysed mesh
% p.Original.Geometrical_features.Volume
% p.Original.Geometrical_features.Principal_inertia_tensor
% 
% % Print some morphological features of the analysed mesh
% p.Original.Morphological_features.Form.Convexity
% p.Original.Morphological_features.Form.Sphericity_Wadell
% p.Original.Morphological_features.Form.Form_indices_obb
% 
% disp('Done!')
