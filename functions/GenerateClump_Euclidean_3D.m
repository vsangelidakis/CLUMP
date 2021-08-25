function [mesh, clump]=GenerateClump_Euclidean_3D( inputGeom, N, rMin, div, overlap, varargin )
%% Clump generator using the Euclidean map for voxelated, 3D particles
% 2021 © V. Angelidakis, S. Nadimi, M. Otsubo, S. Utili.

%% The main concept of this methodology:
% 1. We import the geometry of a particle either as a surface mesh or a
%	 voxelated 3D image.
% 2. If a mesh is given, we transform it into a voxelated representation,
%	 i.e. a binary 3D image, where each voxel belonging to the particle is
%	 equal to zero.
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
%	-inputGeom:	Input geometry, given in one of the formats below:
%				1. Directory of .stl file (for surface meshes)
%				2. Directory of .mat file (for binary voxelated images)
%				3. Struct with fields {vertices,faces} (for surface meshes)
%				4. Struct with fields {img,voxel_size} (for binary voxelated images)
%				   where 
%					- img:			[Nx x Ny x Nz] voxelated image
%					- voxel_size:	[1x3] voxel size in Cartesian space
%
%	- N:		Number of spheres to be generated.
%
%	- rMin:		Minimum allowed radius: When this radius is met, the
%				generation procedure stops even before N spheres are
%				generated.
%
%	- div:		Division number along the shortest edge of the axes-aligned
%				bounding box (AABB) of the particle during voxelisation,
%				i.e. during the transformation of an input surface mesh
%				into a 3D image. It controls the resolution of the
%				voxelated representation of the particle. Not used when a
%				3D image is provided directly as input. If not given,
%				div=50 (default value in iso2mesh).
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
%	- mesh	:	structure containing all relevant parameters of input polyhedron
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
% inputGeom='Hexahedron_Coarse_Mesh.stl'; N=24; rMin=0; div=102; overlap=0.6; output='EU_octaCoarse.txt'; visualise=true;
% [mesh, clump]=GenerateClump_Euclidean_3D( inputGeom, N, rMin, div, overlap, output, visualise );

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
%% Import dependencies
addpath(genpath('../lib'))  % Add path to dependencies (external codes)

%% Configure input particle geometry based on the variable type of inputGeom
switch class(inputGeom)
	case 'char'
		if strcmp(inputGeom(end-3:end),'.stl') % input file is .stl (surface mesh)
			[P,F,~] = stlRead(inputGeom);
			
			% Calculate Rigid Body Parameters (RBP)
			FV=struct();	FV.vertices=P;	FV.faces=F; [RBP,~]=RigidBodyParams(FV);
			
			% Transform surface mesh to voxelated image
			[imgTemp, map]=s2v(P,F,div);
			
			imgTemp2=fillholes3d(imgTemp,2); % This causes some loss of accuracy around the value of 2 voxels (uses imclose), but is needed to ensure that imfill below works properly.
			imgTemp2=imfill(imgTemp2); % Fill the interior of the particle with true values (1).
			% imgTemp2=imfill(imgTemp); % Fill the interior of the particle with true values (1).
			img=zeros(size(imgTemp2)+2); % Expand the image by 2 voxels in each direction, to ensure the boundary voxels are false (zeros).
			img(2:end-1,2:end-1,2:end-1)=imgTemp2;
			clear imgTemp imgTemp2
		
			% Ensure the voxel size is the same in all 3 directions -> Might be an overkill, but still
			if abs((map(1,1)-map(2,2))/map(1,1))<1e-6 || abs((map(2,2)-map(3,3))/map(2,2))<1e-6
				voxel_size=map(1,1);
			else
				warning('The affine transformation from voxels to Cartesian dimensions is not the same in all directions. Voxel size is not the same in X,Y,Z! Potentially wrong radii in Cartesian units!')
				voxel_size=map(1,1);
			end
			
		elseif strcmp(inputGeom(end-3:end),'.mat')  % input file is .mat (voxelated image)
			vox=load(inputGeom);
			temp=fieldnames(vox); temp=temp{1}; vox=vox.(temp); clear temp;
			img=vox.img;
			voxel_size=vox.voxel_size;
			
			opt=2; %see vol2mesh function in iso2mesh
			isovalues=[]; %see vol2mesh function in iso2mesh
			[P,F]=v2s(vox.img,isovalues,opt,'cgalmesh');
			P=P*voxel_size(1,1);
			
			% Calculate Rigid Body Parameters (RBP)
			FV=struct();	FV.vertices=P;	FV.faces=F; [RBP,~]=RigidBodyParams(FV);
		else
			error('Not recognised inputGeom format.')
		end
	case 'struct'
		if isfield(inputGeom,'Vertices') || isfield(inputGeom,'vertices')  % input file is struct containing surface mesh
			try P=inputGeom.Vertices; catch, P=inputGeom.vertices; end
			try F=inputGeom.Faces;    catch, F=inputGeom.faces;	  end

			% Calculate Rigid Body Parameters (RBP)
			[RBP,~]=RigidBodyParams(inputGeom);

			% Transform surface mesh to voxelated image
			[imgTemp, map]=s2v(P,F,div);
			
			imgTemp2=fillholes3d(imgTemp,2); % This causes some loss of accuracy around the value of 2 voxels (uses imclose), but is needed to ensure that imfill below works properly.
			imgTemp2=imfill(imgTemp2); % Fill the interior of the particle with true values (1).
			% imgTemp2=imfill(imgTemp); % Fill the interior of the particle with true values (1).
			img=zeros(size(imgTemp2)+2); % Expand the image by 2 voxels in each direction, to ensure the boundary voxels are false (zeros).
			img(2:end-1,2:end-1,2:end-1)=imgTemp2;
			clear imgTemp imgTemp2			

			% Ensure the voxel size is the same in all 3 directions -> Might be an overkill, but still
			if abs((map(1,1)-map(2,2))/map(1,1))<1e-6 || abs((map(2,2)-map(3,3))/map(2,2))<1e-6
				voxel_size=map(1,1);
			else
				warning('The affine transformation from voxels to Cartesian dimensions is not the same in all directions. Voxel size is not the same in X,Y,Z! Potentially wrong radii in Cartesian units!')
				voxel_size=map(1,1);
			end
			
		elseif isfield(inputGeom,'img')  % input file is struct containing voxelated image
			img=inputGeom.img;
			voxel_size=inputGeom.voxel_size;
			
			opt=2; %see vol2mesh function in iso2mesh
			isovalues=[]; %see vol2mesh function in iso2mesh
			[P,F]=v2s(inputGeom.img,isovalues,opt,'cgalmesh');
			P=P*voxel_size(1,1);

			% Calculate Rigid Body Parameters (RBP)
			FV=struct();	FV.vertices=P;	FV.faces=F; [RBP,~]=RigidBodyParams(FV);
		else
			error('Not recognised inputGeom format.')
		end
	otherwise
		error('Not recognised inputGeom format.')
end

% figure()
% volshow(img)

	%% TODO: Define clump, mesh as classes, not structs.

% The if statement below fixes a bug in the sign of volume/inertia in the RBP code, if the surface normal vectors point inwards.
if RBP.volume<eps
	RBP.volume=-RBP.volume;
% 	RBP.inertia=-RBP.inertia;
% 	RBP.orientationsPrincipal=-1*RBP.orientationsPrincipal;
	disp('Correcting the sign of volume, attributed to inverted normals.') % volume and inertia
end
	
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

%% Dimensions of the new image
halfSize=[size(img,2)/2, size(img,1)/2, size(img,3)/2];

[dx,dy,dz] = meshgrid(1:size(img,2), 1:size(img,1), 1:size(img,3));

%% Calculate centroid of the voxelated image
stats = regionprops3(img,'Centroid'); % 'all'
centroid=stats.Centroid; % Centroid of the initial particle

counter=1;
intersection=img;

for k=1:N %N:numberofspheres
	edtImage = bwdist(~intersection);	% Euclidean map
	radius = max(edtImage(:));			% Inradius in voxel units
	
	% Note: rMin is given in Cartesian units, not in voxel units, hence the
	% multiplication with the voxel size below
	if radius*voxel_size(1,1)<rMin % Break the loop if the minimum radius has been met using less than N spheres
		warning(['The mimimum radius rMin=',num2str(rMin),' has been met using ', num2str(k-1),' spheres'])
		break
	end
	
	[yCenter, xCenter, zCenter]= ind2sub(size(intersection),find(edtImage == radius)); % Center in voxel units
	
	% Sometimes, more than 1 sphere with the same inradius can emerge.
	% Choose the most/least distant sphere to the centroid of the particle.
	%% TODO: Decide whether to use the most/least distant sphere
	%% TODO: In the future, additionally to using the most/least distant sphere, I can use something like Greedy heuristics, to find the sphere with the maximum coverage.
	dists = sqrt(sum(bsxfun(@minus,centroid,[xCenter,yCenter,zCenter]).^2,2));
	[~,i]=max(dists); % Index of the inscribed sphere closest (min) / farthest (max) to the centroid
	
	sph=sqrt( (dx-xCenter(i)).^2 + (dy-yCenter(i)).^2 + (dz-zCenter(i)).^2 ) > (1-overlap)*radius; % Sphere
	intersection=and(intersection,sph); % Append the new sphere in the particle
	
	%% FIXME: Small detail: In some cases we need to add one voxel or subtract one voxel. We need to find a robust criterion on this
	%% FIXME: Small detail: In the beginning, for surface mesh inputs, we enlarged the image by 2 voxels upon its transformation to a 3D voxelated matrix. % When we return to Cartesian space, we need to subtract the length of one voxel's side.
	
	xC=xCenter(i)-halfSize(1); %+1 % Maybe add 1 voxel since we have enlarged the image by 2 for mesh inputs?
	yC=yCenter(i)-halfSize(2); %+1
	zC=zCenter(i)-halfSize(3); %+1
	
	clump.positions(counter,:)=[yC,xC,zC]*voxel_size(1,1)+[aveX,aveY,aveZ]; % Here we add [aveX,aveY,aveZ] to return to the initial coordinate system
	clump.radii(counter,1)=radius*voxel_size(1,1);
	counter=counter+1;
end

[clump.minSphere.centroid, clump.minSphere.radius]=min(clump.radii);
[clump.maxSphere.centroid, clump.maxSphere.radius]=max(clump.radii);
clump.numSpheres=length(clump.radii);

%% Plot spheres in voxelised space (if overlap>0, these are not the actual spheres, but the scaled ones, used to facilitate the overlap)
% figure()
% load('config.mat') % "config" is only used to visualise transparent voxelated images
% volshow(intersection,config);

% Restore the mesh in the original coordinate system
P(:,1)=P(:,1)+aveX;
P(:,2)=P(:,2)+aveY;
P(:,3)=P(:,3)+aveZ;

%% Plot clump and mesh in Cartesian space (optional)
if visualise
% 	patch('Faces',F,'Vertices',P,'FaceColor','g','EdgeColor','none','FaceAlpha',0.4)
% 	patch('Faces',F,'Vertices',P,'FaceColor','g','FaceAlpha',0.15,'EdgeColor',[0.2,0.7,0.2],'EdgeAlpha',0.1);
	patch('Faces',F,'Vertices',P,'FaceColor','g','FaceAlpha',0.2,'EdgeColor','none','EdgeAlpha',0.4); %[0,0.4,0]
	axis equal
	camlight
	box on; grid on
	alpha 0.5
	
	%% Plot spheres
	[X,Y,Z]=sphere(20);
	for i=1:length(clump.radii)
		hold on
		x=clump.positions(i,1);
		y=clump.positions(i,2);
		z=clump.positions(i,3);
		r=clump.radii(i);
		surf(r*X+x, r*Y+y, r*Z+z, 'FaceColor',rand(1,3), 'EdgeColor','none','FaceAlpha',1)
	end
end

%% Export clump (optional)
% Output is offered in the generic format x_y_z_r. 
% For more specialised formats, try the ExportClump module.
if ~isempty(output)
	dlmwrite(output, [clump.positions, clump.radii], 'delimiter', ',', 'precision', '%10f')
end
end
