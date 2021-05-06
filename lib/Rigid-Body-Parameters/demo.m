clc; clear; close all

%% Calculate RBP without any checks
% load('Auxiliary Functions/sample_mesh')
% RBP=RigidBodyParams(TR);
% disp(RBP)
% VisualizeLocalFrame(TR)


%% Calculate RBP checking for the orientation of normals
load('Auxiliary Functions/sample_mesh')
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