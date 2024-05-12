function [mesh, clump]=GenerateClump_Favier( inputGeom, N, varargin )
%% Implementation of the clump-generation concept proposed by Favier et al. (1999) [1]
% 2021 © V. Angelidakis, S. Nadimi, M. Otsubo, S. Utili.
% [1] Favier, J.F., Fard, M.H., Kremmer, M. and Raji, A.O., 1999. Engineering Computations: Int J for Computer-Aided Engineering, 16(4), pp.467-480.

%% The main concept of this methodology:

% 1. We import the geometry of a particle either as a surface mesh or a
%	 voxelated 3D image.
% 2. If a voxelated image is given, we transform it into a surface mesh,
%	 providing its vertices and faces (vertex connectivity).
% 3. We calculate the inertial characteristics of the particle and center
%	 it to its centroid and align it to its principal axes.
% 4. We create a number of N points along the longest particle axis and
%	 identify the particle vertices belonging to each of the newly formed
%	 (N+1) spans.
% 5. We generate a sphere centered to each of the N points. The radius of
%	 each sphere is calculated based on the distances of the vertices
%	 within the span of interest. The default behaviour considers the
%	 minimum distance, although an optional parameter (varargin) exists
%	 that takes the values 'min' (default), 'avg' and 'max'.

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
%	- chooseDistance: Preferred method to specify the radii of the
%					  spheres, which can be either the minimum ('min'), the
%					  average ('avg') or the maximum ('max') distance of
%					  the vertices within the span of interest. 
%
%	- output:	File name for output of the clump in .txt form	(optional)*.
%				If not assigned, a .txt output file is not created.
%
%	-visualise: Whether to plot the clump and mesh (optional)*.
%
% * varargin can contain the 'chooseDistance', the 'output' and the
%	'visualise' variables. They are all optional.

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

%% TODO: Maybe it would be helpful to save the centroid/volume/inertia of the whole clump for uniform/non-uniform density
%%
%				clump.centroid		:	Centroid of clump assuming uniform
%										(or non-uniform) density for all spheres
%				clump.volume		:	Volume of clump before (or after)
%										density correction
%				clump.inertia		:	Inertia of clump before (or after)
%										density correction
%				clump.volumes		:	Volume of each sphere
%
%				clump.inertias		:	Inertia of each sphere
%

%% EXAMPLE
% 	inputGeom='ParticleGeometries/Cylinder.stl'; N=20; chooseDistance='min'; output='FA_cylinder.txt'; visualise=true;
% 	[mesh, clump]=GenerateClump_Favier( inputGeom, N, chooseDistance, output, visualise );

%% TODO
% We can calculate the centroid/volume/inertia of the clump and compare it to that of the actual particle
% Ensure the vertex normals are pointing inwards for concave meshes (so far this is done only for convex meshes)

%% Define variables based on the type of the optional parameters (varargin)
chooseDistance='min'; % Default method to choose radius.
for i=1:length(varargin)
	switch class(varargin{i})
		% 		case 'double'
		% 			seed=varargin{i}; rng(seed);	% Fixed seed to achieve reproducible (random) results
		case 'char'
			if strcmp(varargin{i},'min') || strcmp(varargin{i},'avg') || strcmp(varargin{i},'max')
				chooseDistance=varargin{i};
			else
				output=varargin{i};
			end
		case 'logical'
			visualise=varargin{i};
		otherwise
			error('Wrong optional parameter type.')
	end
end

%% Main body of the function
%% Import Dependencies
addpath(genpath('../lib')) % Add path to dependencies (external codes)

%% Configure input particle geometry based on the variable type of inputGeom
switch class(inputGeom)
	case 'char'
		if strcmp(inputGeom(end-3:end),'.stl') % input file is .stl (surface mesh)
			[P,F,~] = stlRead(inputGeom);
		elseif strcmp(inputGeom(end-3:end),'.mat')  % input file is .mat (voxelated image)
			vox=load(inputGeom);
			temp=fieldnames(vox); temp=temp{1}; vox=vox.(temp); clear temp;
			voxel_size=vox.voxel_size;
% 			vox

			% TODO: Add error, if the fields: vox.img and vox.voxel_size do not exist
			
			opt=2; %see vol2mesh function in iso2mesh
			isovalues=[]; %see vol2mesh function in iso2mesh
			[P,F]=v2s(vox.img,isovalues,opt,'cgalmesh');
			P=P*voxel_size(1,1);
		else
			error('Not recognised inputGeom format.')
		end
	case 'struct'
		if isfield(inputGeom,'Vertices') || isfield(inputGeom,'vertices')  % input file is struct containing surface mesh
			try P=inputGeom.Vertices; catch, P=inputGeom.vertices; end
			try F=inputGeom.Faces;    catch, F=inputGeom.faces;	  end
		elseif isfield(inputGeom,'img')  % input file is struct containing voxelated image
			voxel_size=inputGeom.voxel_size;
			opt=2; %see vol2mesh function in iso2mesh
			isovalues=[]; %see vol2mesh function in iso2mesh
			[P,F]=v2s(inputGeom.img,isovalues,opt,'cgalmesh');
			P=P*voxel_size(1,1);
		else
			error('Not recognised inputGeom format.')
		end
	case 'triangulation'
		try
			F=inputGeom.ConnectivityList;
			P=inputGeom.Points;
		catch
			error('Not recognised -triangulation- format.')
		end
	otherwise
		error('Not recognised inputGeom format.')
end

% Ensure all face normals are oriented coherently, pointing outwards
TR2=triangulation(F,P);
[TR,~]=ConsistentNormalOrientation(TR2); %numInvFaces

[RBP,~]=RigidBodyParams(TR);
% Attention: For cubic particles, with no elongation, the RigidBodyParams
% function mis-identifies the principal planes as the ones of the diagonal.
% This is the same mistake the PCA typically does for particles with three
% equal dimensions. Thankfully, the method of Favier et al (1999) is meant
% to be used for elongated and somewhat axi-symmetric particles.

% %	Plot original particle
% 	patch('Faces',F,'Vertices',P,'FaceColor','r') %'r'
% 	axis equal; camlight

%% Center particle around its centroid and align it along its principal axes.
rot=RBP.PAI;

% Transform rotation matrix to align longest axis along X direction
temp=rot(:,1);
rot(:,1)=rot(:,3);
rot(:,3)=temp;

P=P-RBP.centroid;
P=P*rot;

% %	Plot particle centered around its centroid and aligned to its principal axes
% 	patch('Faces',F,'Vertices',P,'FaceColor','g','FaceAlpha',0.3) %'r'
% 	axis equal; camlight

X_extremas=[min(P(:,1)),max(P(:,1))];

a=X_extremas(1);
b=X_extremas(2);
nSegments=N;
%Example:
endPoints = linspace(a,b,nSegments + 1);   %4 endpoints for 3 segments
start = endPoints(1:end-1);                %3 starting points
stop = endPoints(2:end);                   %3 stopping points
midPoints = stop - ((stop(1)-start(1))/2); %3 middle points

% 	scatter(endPoints,zeros(length(endPoints)),'r','filled') % endPoints
% 	scatter(midPoints,zeros(length(midPoints)),'b','filled') % midPoints

minDistance=zeros(1,length(midPoints));
minDx=zeros(1,length(midPoints));
for i=1:length(midPoints) % For each midpoint
	count=1;
	% Find vertices within each sector
	for j=1:size(P,1) % for each point on the particle surface (in principal axes)
		if P(j,1)>=endPoints(i) && P(j,1)<=endPoints(i+1)
			p{i}(count,1:3)=P(j,1:3);
			count=count+1;
		end
	end
	% 	scatter3(p{1,i}(:,1),p{1,i}(:,2),p{1,i}(:,3),20,i/length(midPoints)*rand(1,3),'filled'); hold on %20*i
	
	% Find closest distance of midpoint to any surface vertex
	xM=midPoints(i);   yM=0; zM=0;
	x1=endPoints(1);   %y1=0; z1=0;
	x2=endPoints(end); %y2=0; z2=0;
	
	% Closest distance between midpoint and particle surface
	minDistance(i)=min(sqrt( (P(:,1)-xM).^2 + (P(:,2)-yM).^2 + (P(:,3)-zM).^2 ));
	
	% Closest distane between midpoint and particle X limits
	minDx(i)=min( abs(x1-xM) , abs(x2-xM) );
	
	% 	minDx(i)=min(sqrt( (x1-xM).^2 + (y1-yM).^2 + (z1-zM).^2 ),...
	% 		sqrt( (x2-xM).^2 + (y2-yM).^2 + (z2-zM).^2 ) );
end

% TODO: Check if neighbouring particles include one another (one is
% almost fully inside the other) and if yes, delete the smaller one.
% Then, we would end up with a smaller number of spheres than the
% target one, and so we should run more iterations of the code, using
% a target of (N+1) generation points, until the target number of
% spheres is achieved.

%% Build "clump" structure
clump=struct;
clump.positions=[];
clump.radii=[];

% 	radius=zeros(1,midPoints);
radius=zeros(1,length(midPoints));
for i=1:length(midPoints)
	if ~isempty(p{1,i})
		distance=sqrt( (p{1,i}(:,1)-midPoints(i)).^2 + (p{1,i}(:,2)-0).^2 + (p{1,i}(:,3)-0).^2 );
		switch chooseDistance
			case 'min'
				radius(i)=min(distance);	% Minimum distance
			case 'avg'
				radius(i)=mean(distance);	% Average distance
			case 'max'
				radius(i)=max(distance);	% Maximum distance
 		%	otherwise
 		% 		radius(i)=					% Radius of best fitting sphere TODO
		end
		% [min(distance),mean(distance),max(distance)]
		% Limit radius to not exceed the particle length (along X axis)
		radius(i)=min(radius(i),minDx(i));
		% TODO: Add optional boolean whether to trim the corner spheres or not
		
		% The code below uses the minimum distance to the particle
		% surface, not the smallest to the particle X limits
		% 		radius(i)=min(radius(i),minDistance(i));
		% 		radius(i)=minDistance(i);
		
		clump.positions(i,:)=[midPoints(i),0,0];
		clump.radii(i,1)=radius(i);
	else
		error('The number of particle vertices is small for the requested number of spheres. Either remesh the particle surface or choose a smaller number of spheres.')
		% TODO: Find a more elegant solution if we use so many spheres that we don't have any vertices within the sector we study
	end
end

%% Transform the mesh and the clump coordinates back to the initial (non-principal) system
P=P*rot'; % use either the transposed rot' or inverse inv(rot) rotation matrix to return to the initial coordinate system
P=P+RBP.centroid;

clump.positions=clump.positions*rot';
clump.positions=clump.positions+RBP.centroid;

%% Build "mesh" and finalise "clump" structures
mesh=struct;
mesh.vertices=P;
mesh.faces=F;
mesh.centroid=RBP.centroid;
mesh.volume=RBP.volume;
mesh.inertia=RBP.inertia_tensor;
mesh.inertiaPrincipal=RBP.eigs;
mesh.orientationsPrincipal=RBP.PAI;

% TODO: Define clump, mesh as classes, not structs

[clump.minSphere.centroid, clump.minSphere.radius]=min(clump.radii);
[clump.maxSphere.centroid, clump.maxSphere.radius]=max(clump.radii);
clump.numSpheres=length(clump.radii);

%% Plot clump and mesh (optional)
if visualise
	patch('Faces',F,'Vertices',P,'FaceColor','g','FaceAlpha',0.5,'EdgeColor','none');
	axis equal
	camlight
	box on; grid on; hold on
	alpha 0.4
	
	%% Plot spheres
	for j=1:length(clump.radii)
		[X,Y,Z]=sphere;
		xSph=X*clump.radii(j);
		ySph=Y*clump.radii(j);
		zSph=Z*clump.radii(j);
		
		color=rand(1,3);
		xC=clump.positions(j,1);
		yC=clump.positions(j,2);
		zC=clump.positions(j,3);
		
		surf(xSph+xC,ySph+yC,zSph+zC,'EdgeColor','none','FaceAlpha',0.6,'FaceColor',color)
	end
end

%% Export clump (optional)
% Output is offered in the generic format x_y_z_r. For more specialised
% formats, try the exportClump module.
if isempty(output)==false
	dlmwrite(output, [clump.positions, clump.radii], 'delimiter', ',', 'precision', '%10f')
end

end
