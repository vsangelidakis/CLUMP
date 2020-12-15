%% Clump generator using Euclidean maps for 3D particles
% Copyright © 2020 V. Angelidakis. All rights reserved
tic
clc; clear; close all

addpath(genpath('functions'))

%% INPUT
N=80; % Number of spheres to be generated.
rMin=0; % Minimum allowed radius: Larger values will break the loop even before N spheres are generated.
div=102; % Division number along the shortest edge of the AABB (resolution). If not given, div=50 (default value in iso2mesh).
overlap=0.7; % Overlap percentage: [0,1): 0 for non-overlapping spheres, 0.4 for 40% overlap of radii, etc.


%% Load particle shape from stl
% [P,F,n] = stlRead('Hexahedron_Fine_Mesh.stl');
% [P,F,n] = stlRead('Cone.stl');

% [P,F,n] = stlRead('Cube_3_2_1.stl');
% [P,F,n] = stlRead('Cube_3_2_1_moved.stl');

[P,F,n] = stlRead('stl_stl_5000_faces.stl');
% [P,F,n] = stlRead('LB2_Original_voxel_scale_13078_faces.stl');
% [P,F,n] = stlRead('LB3_Original_voxel_scale_9248_faces.stl');


%% Calculate extreme coordinates & centroid of the AABB of the particle
minX=min(P(:,1)); maxX=max(P(:,1)); aveX=mean([minX,maxX]); %ave: centroid of the AABB
minY=min(P(:,2)); maxY=max(P(:,2)); aveY=mean([minY,maxY]);
minZ=min(P(:,3)); maxZ=max(P(:,3)); aveZ=mean([minZ,maxZ]);

%% Center the particle to the centroid of its AABB
P(:,1)=P(:,1)-aveX;
P(:,2)=P(:,2)-aveY;
P(:,3)=P(:,3)-aveZ;

[imgTemp, map]=s2v(P,F,div);

imgTemp2=fillholes3d(imgTemp,2); % This causes some loss of accuracy around the value of 2 voxels (uses imclose), but is needed to ensure that imfill below works properly.
imgTemp2=imfill(imgTemp2); % Fill the interior of the particle with true values (1).
% imgTemp2=imfill(imgTemp); % Fill the interior of the particle with true values (1).
img=zeros(size(imgTemp2)+2); % Expand the image by 2 voxels in each direction, to ensure the boundary voxels are false (zeros).
img(2:end-1,2:end-1,2:end-1)=imgTemp2;
clear imgTemp imgTemp2

% stats=regionprops3(img,'Volume'); %all

% %% Calculate centroid of the original particle
% TR.vertices=P;
% TR.faces=F;
% 
% [RBP,TR]=RigidBodyParams(TR);
% disp(RBP.centroid)

%% 
if abs((map(1,1)-map(2,2))/map(1,1))<1e-6 || abs((map(2,2)-map(3,3))/map(2,2))<1e-6
	voxel_size=map(1,1);
else
	warning('The affine transformation from voxels to Cartesian dimensions is not the same in all directions. Voxel size is not the same in X,Y,Z! Potentially wrong radii in Cartesian units!')
	voxel_size=map(1,1);
end

%% Alternatively, uncomment to create cuboidal particle directly in voxelised space
% img=false(102,102,202);		% Create domain
% img(10:90, 10:90, 10:190)=true;	% Create cuboidal particle
% voxel_size=0.02;
%% Hollow clumps, e.g. shelly carbonate sand?
% img(35:65, 35:65, 75:125)=false;	% Create an inner layer for hollow clumps
% Hollow clumps can be achieved by setting the interior of the particle as
% false. The user has to decide how hollow, i.e. what is the diameter of
% the main spheres to be used.
%%

%% Dimensions of the new image
halfSize=[size(img,2)/2, size(img,1)/2, size(img,3)/2];

[dx,dy,dz] = meshgrid(1:size(img,2), 1:size(img,1), 1:size(img,3));

%% Calculate centroid of the voxelated image
stats = regionprops3(img,'Centroid'); % 'all'
centroid=stats.Centroid; % Centroid of the initial particle

counter=1;
ind=1;
spheresList=[];
intersection=img;

for k=1:N
	%%
	edtImage = bwdist(~intersection);		% Euclidean map
	radius = max(edtImage(:));				% Inradius in voxel units
	
	% Break the loop if the minimum radius has been met using less than N spheres
	if radius<rMin
		disp(['The mimimum radius rMin=',num2str(rMin),' has been met using ', num2str(k-1),' spheres'])
		break
	end
	
	[yCenter, xCenter, zCenter]= ind2sub(size(intersection),find(edtImage == radius)); % Center in voxel units
	
	% Sometimes, more than 1 sphere with the same inradius can emerge
	% Find the most/least distant sphere to the centroid of the particle
	% FIXME: Need to decide whether to use the most/least distant sphere
	dists = sqrt(sum(bsxfun(@minus,centroid,[xCenter,yCenter,zCenter]).^2,2));
	[~,i]=max(dists); % Index of the inscribed sphere closer (min) / farther (max) to the centroid
	
	sph=zeros(length(dy),length(dx),length(dz));
	sph=sqrt( (dx-xCenter(i)).^2 + (dy-yCenter(i)).^2 + (dz-zCenter(i)).^2 ) > (1-overlap)*radius; % Sphere
	intersection=and(intersection,sph);

	%% FIXME: I need to remove the spheres that are fully inside the particle and thus play no role. 
	%% I can check whether the extreme voxels of every new sphere are touching the 
	%% boundaries of the initial particle or alternatively, if the extreme voxels 
	%% of each new sphere are farther than the extreme voxels of the existing spheres.

	%% FIXME: In some cases we need to add one voxel or subtract one voxel. I need to find a robust criterion on this
	%% In the beginning, I enlarged the image by 2 voxels. % When I return to Cartesian space, I need to subtract the length of one voxel's side.
	
	xC=xCenter(i)-halfSize(1)+1; % I add 1 voxel since we have enlarged the image by 1
	yC=yCenter(i)-halfSize(2)+1;
	zC=zCenter(i)-halfSize(3)+1;

% 	spheresList(ind,1:4)=[yCenter(i),-xCenter(i),zCenter(i),radius]*voxel_size; % Cartesian units
	spheresList(ind,1:4)=[yC,xC,zC,radius]*voxel_size; % Cartesian units
	ind=ind+1;

	
end

%% FIXME: I need to bring the actual particle and the clumps back to their initial coordinate system, using aveX,Y,Z!

% 	newno=map*[no ones(length(no),1)]';
% 	newno=newno(1:3,:)'; % newno and el now go back to the world coordinates

%% Plot spheres in voxelised space (if overlap>0, these are not the actual spheres, but the scaled ones, used to facilitate the overlap)
figure()
load('config.mat')
volshow(intersection,config);

%% Plot clump in Cartesian space
[X,Y,Z]=sphere(40);
figure()
for i=1:size(spheresList,1)
	hold on
	x=spheresList(i,1);
	y=spheresList(i,2);
	z=spheresList(i,3);
	r=spheresList(i,4);
	surf(r*X+x, r*Y+y, r*Z+z, 'FaceColor',rand(1,3), 'EdgeColor','none','FaceAlpha',0.7)
end
axis equal
grid on
camlight

patch('Faces',F,'Vertices',P,'FaceColor','g','EdgeColor','none','FaceAlpha',0.4)

disp([newline,'Done!'])
toc
