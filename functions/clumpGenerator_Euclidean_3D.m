function [mesh, clump]=clumpGenerator_Euclidean_3D( stlFile, N, rMin, div, overlap, varargin )

%% Clump generator using the Euclidean map for voxelated, 3D particles
% 2020 © V. Angelidakis, S. Nadimi, M. Otsubo, S. Utili.

%% The main concept of this methodology:
% 1. We import the surface mesh of a particle.
% 2. We transform the mesh into a voxelated representation, i.e. a binary
%    3D image, where each voxel belonging to the particle is equal to zero.
% 3. The Euclidean distance transform of the 3D image is computed and
%    the radius of the largest inscribed sphere is found as the maximum
%    value of the Euclidean transform of the voxelated image.
% 4. The voxels corresponding to the inscribed sphere are then set equal to
%    one. This methodology can also generate overlapping spheres, if only a
%    percentage of the voxels of each new sphere are set equal to one,
%    instead of all of them.
% 5. This process is repeated until a user-defined number of spheres 'N' is
%    found or until the user-defined minimum radius criterion has been met,
%    as the spheres are generated in decreasing sizes.

%% Influence of parameters
% N:		[1,inf)  Larger N will lead to a larger number of spheres
% rMin:		(0,inf)  Larger rMin will lead to a smaller number of spheres
% div:		(5,inf]  Larger div will lead to better shape resolution in voxel space
% overlap:	[0,1)	 Larger overlap will lead to larger spheres overall

%% INPUT:
%	- stlFile:	Directory of stl file, used to generate spheres
%
%	- N:		Number of spheres to be generated.
%
%	- rMin:		Minimum allowed radius: When this radius is met, the
%				generation procedure stops even before N spheres are
%				generated.
%
%	- div:		Division number along the shortest edge of the AABB during
%				voxelisation (resolution). If not given, div=50 (default
%				value in iso2mesh).
%
%	- overlap:	Overlap percentage: [0,1): 0 for non-overlapping spheres,
%				0.4 for 40% overlap of radii, etc.
%
%	- output:	File name for output of the clump in .txt form (optional)*.
%				If not assigned, a .txt output file is not created.
%
%	-visualise: Whether to plot the clump and mesh (optional)*.
%
% * varargin can contain either of the optional variables "output",
% "visualise" or else: output=varargin{1}; visualise=varargin{2}.

%% OUTPUT:
%	- mesh	:	structure containing all relevant parameters of polyhedron
%				mesh.vertices
%				mesh.faces
%				mesh.centroid
%				mesh.volume
%				mesh.inertia
%				mesh.inertiaPrincipal
%				mesh.orientationsPrincipal
%
%   - clump	:	structure containing all relevant clump parameters
%               clump.positions		:	M-by-3 matrix containing the
%										position of each generated sphere.
%				clump.radii			:	M-by-1 vector containing the radius
%										of each generated sphere
%				clump.minRadius		:	Minimum generated sphere (might
%										differ from rmin)
%				clump.maxRadius		:	Maximum generated sphere
%
%				clump.numSpheres	:	Total number of spheres
%
%	- output :	txt file with centroids and radii, with format: [x,y,z,r]

%% EXAMPLE
% stlFile='Hexahedron_Fine_Mesh.stl'; N=24; rMin=0; div=102; overlap=0.6; output='EU_octaCoarse.txt'; visualise=true;
% [mesh, clump]=clumpGenerator_Euclidean_3D( stlFile, N, rMin, div, overlap, output, visualise );

%% Define variables based on the type of the optional parameters (varargin)
output=[];
visualise=false;
for i=1:length(varargin)
	switch class(varargin{i})
		case 'char'
			output=varargin{i};
		case 'logical'
			visualise=varargin{i};
		otherwise
			error('Wrong optional parameter type.')
	end
end

%% Main body of the function
%% Import Dependencies
addpath(genpath('../lib'))  % Add path to dependencies (external codes)

[P,F,~] = stlRead(stlFile);

%% FV structure with faces/vertices, to be used in RBP below
FV=struct();
FV.vertices=P;
FV.faces=F;

[RBP,~]=RigidBodyParams(FV);

%% TODO: Define clump, mesh as classes, not structs.

%% Build "mesh" structure
mesh=struct;
mesh.vertices=P;
mesh.faces=F;
mesh.centroid=RBP.centroid;
mesh.volume=RBP.volume;
mesh.inertia=RBP.inertia_tensor;
mesh.inertiaPrincipal=RBP.eigs;
mesh.orientationsPrincipal=RBP.PAI;

%% Build "clump" structure
clump=struct;
clump.positions=[];
clump.radii=[];

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

% %% Calculate centroid of the original particle
% TR.vertices=P;
% TR.faces=F;
%
% [RBP,TR]=RigidBodyParams(TR);
% disp(RBP.centroid)

%% Ensure the voxel size is the same in all 3 directions -> Might be an overkill, but still
if abs((map(1,1)-map(2,2))/map(1,1))<1e-6 || abs((map(2,2)-map(3,3))/map(2,2))<1e-6
	voxel_size=map(1,1);
else
	warning('The affine transformation from voxels to Cartesian dimensions is not the same in all directions. Voxel size is not the same in X,Y,Z! Potentially wrong radii in Cartesian units!')
	voxel_size=map(1,1);
end


%% Dimensions of the new image
halfSize=[size(img,2)/2, size(img,1)/2, size(img,3)/2];

[dx,dy,dz] = meshgrid(1:size(img,2), 1:size(img,1), 1:size(img,3));

%% Calculate centroid of the voxelated image
stats = regionprops3(img,'Centroid'); % 'all'
centroid=stats.Centroid; % Centroid of the initial particle

counter=1;
intersection=img;

for k=1:N
	edtImage = bwdist(~intersection);		% Euclidean map
	radius = max(edtImage(:));				% Inradius in voxel units
	
	if radius<rMin % Break the loop if the minimum radius has been met using less than N spheres
		warning(['The mimimum radius rMin=',num2str(rMin),' has been met using ', num2str(k-1),' spheres'])
		break
	end
	
	[yCenter, xCenter, zCenter]= ind2sub(size(intersection),find(edtImage == radius)); % Center in voxel units
	
	% Sometimes, more than 1 sphere with the same inradius can emerge.
	% Choose the most/least distant sphere to the centroid of the particle.
	%% FIXME: Decide whether to use the most/least distant sphere
	%% TODO: In the future, additionally to using the most/least distant sphere, I can use something like Greedy heuristics, to find the sphere with the maximum coverage.
	dists = sqrt(sum(bsxfun(@minus,centroid,[xCenter,yCenter,zCenter]).^2,2));
	[~,i]=max(dists); % Index of the inscribed sphere closer (min) / farther (max) to the centroid
	
	sph=sqrt( (dx-xCenter(i)).^2 + (dy-yCenter(i)).^2 + (dz-zCenter(i)).^2 ) > (1-overlap)*radius; % Sphere
	intersection=and(intersection,sph);
	
	%% FIXME: In some cases we need to add one voxel or subtract one voxel. I need to find a robust criterion on this
	%% In the beginning, I enlarged the image by 2 voxels. % When I return to Cartesian space, I need to subtract the length of one voxel's side.
	
	xC=xCenter(i)-halfSize(1)+1; % I add 1 voxel since we have enlarged the image by 1
	yC=yCenter(i)-halfSize(2)+1;
	zC=zCenter(i)-halfSize(3)+1;
	
	clump.positions(counter,:)=[yC,xC,zC]*voxel_size;
	clump.radii(counter,1)=radius*voxel_size;
	counter=counter+1;
end

%% FIXME: Test if I need this!
% clump.positions(:,1)=clump.positions(:,1)+aveX;
% clump.positions(:,2)=clump.positions(:,2)+aveY;
% clump.positions(:,3)=clump.positions(:,3)+aveZ;

[clump.minSphere.centroid, clump.minSphere.radius]=min(clump.radii);
[clump.maxSphere.centroid, clump.maxSphere.radius]=max(clump.radii);
clump.numSpheres=length(clump.radii);

%% FIXME: I need to bring the actual particle and the clumps back to their initial coordinate system, using aveX,Y,Z!
%% Potentially add a bool, whether to center the "mesh" and "clump" to the centroid of the mesh or else use the initial, random, Coordinate System

% 	newno=map*[no ones(length(no),1)]';
% 	newno=newno(1:3,:)'; % newno and el now go back to the world coordinates

%% Plot spheres in voxelised space (if overlap>0, these are not the actual spheres, but the scaled ones, used to facilitate the overlap)
% figure()
% load('config.mat') % "config" is only used to visualise transparent voxelated images
% volshow(intersection,config);

%% Plot clump and mesh in Cartesian space
if visualise
	[X,Y,Z]=sphere(40);
	% figure()
	for i=1:length(clump.radii)
		hold on
		x=clump.positions(i,1);
		y=clump.positions(i,2);
		z=clump.positions(i,3);
		r=clump.radii(i);
		surf(r*X+x, r*Y+y, r*Z+z, 'FaceColor',rand(1,3), 'EdgeColor','none','FaceAlpha',0.7)
	end
	axis equal
	grid on
	camlight
	
	patch('Faces',F,'Vertices',P,'FaceColor','g','EdgeColor','none','FaceAlpha',0.4)
end

%% Export clump
% Output is offered in the generic format x_y_z_r. For more specialised
% formats, try the exportClump module.
if isempty(output)==false
	dlmwrite(output, [clump.positions, clump.radii], 'delimiter', ',', 'precision', '%10f')
end
end
