function [mesh, clump]=clumpGenerator_Ferellec_McDowell( stlFile, dmin, rmin, rstep, pmax, varargin )
%% Implementation of the clump-generation concept proposed by Ferellec and McDowell (2010) [1]
% 2020 © V. Angelidakis, S. Nadimi, M. Otsubo, S. Utili.
% [1] Ferellec, J.F. and McDowell, G.R., 2010. Granular Matter, 12(5), pp.459-467. DOI 10.1007/s10035-010-0205-8

%% The main concept of this methodology:
% 1. We import the surface mesh of a particle.
% 2. We calculate the normal of each vertex pointing inwards.
% 3. For a random vertex on the particle surface, we start creating
%	 tangent spheres with incremental radii along the vertex normal,
%	 starting from 'rmin', with a step of 'rstep', until they meet the
%	 surface of the particle.
% 4. We select a new vertex randomly, which has a distance larger
%    than 'dmin' from the existing spheres and do the same.
% 5. When a percentage 'pmax' of all the vertices is used to generate
%    spheres, the generation procedure stops.
% -  An optional 'seed' parameter is introduced, to generate reproducible
%    clumps.

%% Influence of parameters
% rmin:	(0,inf) Larger rmin will lead to a smaller number of spheres
% dmin: [0,inf) Larger dmin will lead to a smaller number of spheres
% pmax: (0,1]   Larger pmax will lead to a larger number of spheres

% Pros:		The authors of this methodology claim efficiency and preservation
%			of flat faces (reduced artificial roughness compared to other techniques).
% Cons:		The methodology is mesh-dependent, as spheres are generated at the
%			vertices of the input mesh.
% Warning:	The authors of this methodology advise that if the initial mesh
%			is very finely discretised, an adequately large rmin value
%			should be used, to guard the process against "parasitic
%			spheres", i.e. very small spheres which might result to
%			numerical instabilities when using DEM.

%% INPUT:
%	- stlFile:	Directory of stl file, used to generate spheres
%
%	- dmin	:	Minimum allowed distance between new vertex of the surface
%				mesh and existing spheres. If left zero, this distance is
%				not cheched.
%
%	- rmin	:	Minimum radius of sphere to be generated. For coarse
%				meshes, the actual minimum radius might be >rmin.
%
%	- rstep	:	Step used to increase the radius in each iteration, until
%				the generated sphere meets another point of the particle.
%
%	- pmax	:	Percentage of vertices which will be used to generate
%				spheres. The selection of vertices is random.
%
%	- seed	:	Seed value, used to achieve reproducible (random) results (optional)*
%
%	- output:	File name for output of the clump in .txt form	(optional)*.
%				If not assigned, a .txt output file is not created.
%
%	-visualise: Whether to plot the clump and mesh (optional)*.
%
% * varargin can contain either of the optimal variables "seed", "output",
% visualise" or else: seed=varargin{1}; output=varargin{2}; visualise=varargin{3}.

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
% stlFile='Hexahedron_Fine_Mesh.stl';	dmin=0.01;	rmin=0.01; rstep=0.001;	pmax=1.0;	seed=5;	output='hexaFine.txt'; visualise=true;
% [mesh, clump]=clumpGenerator_Ferellec_McDowell( stlFile, dmin, rmin, rstep, pmax, seed, output, visualise );

%% TODO
% We can calculate the centroid/volume/inertia of the clump and compare it to that of the actual particle
% Ensure the vertex normals are pointing inwards for concave meshes (so far this is done only for convex meshes)

%% Define variables based on the type of the optional parameters (varargin)
output=[];
visualise=false;
for i=1:length(varargin)
	switch class(varargin{i})
		case 'double'
			seed=varargin{i}; rng(seed);	% Fixed seed to achieve reproducible (random) results
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
addpath(genpath('../lib')) % Add path to dependencies (external codes)

[P,F,~] = stlRead(stlFile);

FV=struct();
FV.vertices=P;
FV.faces=F;

N=patchnormals(FV);	%% FIXME: Revisit to replace this with Matlab's in-house function to calculate normals of triangulations. To this, I can transform FV into a triangulation object, instead of a "patch" style struct. Also, remove patchnormals from the dependencies.
[RBP,~]=RigidBodyParams(FV);
% VisualizeLocalFrame(TR,3)

%% Ensure the vertex normals are pointing inwards
%% ATTENTION: This approach is guaranteed to work only for convex polyhedra with manifold meshes.
% For concave polyhedra, the vector of a random vertex to the centroid
% does not necessarily reflect whether a face is poining inwards or outwards.

% For concave particles, we can generate a tetrahedral mesh of the
% particle from the surface mesh and use one of the non-parallel edges
% of the adjacent tetrahedron to specify the inwards direction.

for i=1:size(P,1)
	if dot(P(i,:)-RBP.centroid,N(i,:))>0
		N(i,:)=-N(i,:);
	end
end

Pmax=1:length(P);	% List of vertices indices

% if nargin>5
% 	if ischar(varargin{1})==false && isempty(varargin{1})==false
% 		rng(seed)		% Fixed seed to achieve reproducible (random) results
% 	end
% end

Vertices=Pmax(randperm(length(Pmax)));	% Shuffle indices of vertices (random selection)
% Vertices=Pmax;						% Ordered indices of vertices (ordered selection)

tol=rmin/1000;	% Tolerance so that the starting vertex is considered outside the sphere


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

counter=1;
iCount=1;
for k=Pmax
	i=Vertices(iCount);
	r=rmin;
	reachedMaxRadius=false;	% The maximum radius is reached when the sphere becomes large enough to include a point of the mesh
	
	x=P(i,1); % Vertex used to generate sphere
	y=P(i,2);
	z=P(i,3);
	
	n=N(i,:); % Vertex normal facing inwards
	
	%% Check if vertex is closer than dmin to the surface of one of the existing spheres
	if iCount>1 && dmin>0
		dcur=min(sqrt( (x-clump.positions(:,1)).^2 + (y-clump.positions(:,2)).^2 + (z-clump.positions(:,3)).^2 ) - clump.radii(:) );	% Distances of all points to the center of the sphere
		if dcur<dmin
			iCount=iCount+1;
			continue
		end
	end
	
	%% Alternative, non-vectorised loop for dmin check
	% 		skipVertex=false;
	% 		for m=1:length(clump.radii)
	% 			dcur=min(sqrt( (x-clump.positions(m,1)).^2 + (y-clump.positions(m,2)).^2 + (z-clump.positions(m,3)).^2 ) - clump.radii(m) );
	% 			if abs(dcur)<clump.radii(m) && dcur<dmin
	% 				skipVertex=true;
	% 			end
	% 		end
	%
	% 		if skipVertex
	% 			iCount=iCount+1;
	% 			continue
	% 		end
	
	% 	scatter3(x, y, z,'r') % Uncomment to visualise each point that is used to generate spheres
	
	while reachedMaxRadius==false % while the sphere has not reached the particle surface
		sphMin=1e15; % Minimum value of potential function
		while sphMin>-tol
			xC=x+r*n(1);
			yC=y+r*n(2);
			zC=z+r*n(3);
			
			distance=sqrt( (P(:,1)-xC).^2 + (P(:,2)-yC).^2 + (P(:,3)-zC).^2 );	% Distances of all points to the center of the sphere
			sph=(distance/r).^2-1;												% Value of spherical potential function (negative for points inside the sphere)
			sphMin=min(sph);
			
			r=r+rstep; % Grow radius for next step
			
			%% FIXME: Maybe instead of calculating the distance to all vertices, calculate the distance to all faces!!!
			
			%% TODO: Maybe here check that if the normal is pointing outwards, we inverse the sign of n(1), n(2), n(3).
			%% This can be done either when we reach a maximum number of iterations or more safely, if the distance decreases instead of increasing, for every iteration.
		end
		reachedMaxRadius=true;
		indMin=find(sph==sphMin);
		
		pointInside=P(indMin(1),:);
		
		% 		for k=1:length(indMin) % Uncomment to visualise the points of the mesh that are included in the current sphere
		% 			scatter3(P(indMin(k),1), P(indMin(k),2), P(indMin(k),3),'b')
		% 		end
		
		vAB=[pointInside(1)-x, pointInside(2)-y, pointInside(3)-z]; % Vector from starting point to point with min distance to the center of the sphere
		vAD=dot(vAB,n)/norm(n);										% Projection of previous vector on the current normal vector
		% 		theta =atan2( vecnorm(cross(n,vAB,2),2,2) , dot(n,vAB,2) ); % Angle between vAB and normal vector n
		
		AB=norm(vAB);
		AD=norm(vAD);
		% 		BD=sqrt(AB^2-AD^2);
		
		radius=AB^2/AD/2;
		
		xC=x+radius*n(1);
		yC=y+radius*n(2);
		zC=z+radius*n(3);
		
		clump.positions(counter,:)=[xC,yC,zC];
		clump.radii(counter,1)=radius;
		counter=counter+1;
		
	end
	%% Check whether the maximum percentage of vertices has been used
	pcur=length(clump.radii)/length(P); % Current percentage of vertices used
	if pcur<pmax
		iCount=iCount+1;
	else
		break
	end
end

[clump.minSphere.centroid, clump.minSphere.radius]=min(clump.radii);
[clump.maxSphere.centroid, clump.maxSphere.radius]=max(clump.radii);
clump.numSpheres=length(clump.radii);

%% Plot clump and mesh
if visualise
	patch('Faces',F,'Vertices',P,'FaceColor','g','FaceAlpha',0.5,'EdgeColor','none') ;%[0.5,0.5,0.5]
	axis equal
	camlight
	% hL1=camlight('headlight');
	% set(hL1,'style','infinite','position',mesh.centroid*2)
	% set(gca,'visible','off')
	box on; grid on; hold on
	alpha 0.5
	
	%% Plot normals (they should point inwards)
	% quiver3(P(:,1),P(:,2),P(:,3),N(:,1),N(:,2),N(:,3)) % Uncomment to visualise normal vectors of vertices. They should all point inwards
	
	%% Plot spheres
	for j=1:length(clump.radii)
		[X,Y,Z]=sphere;
		xSph=X*clump.radii(j);
		ySph=Y*clump.radii(j);
		zSph=Z*clump.radii(j);
		
		color=rand(1,3);
		% 	color='g';
		xC=clump.positions(j,1);
		yC=clump.positions(j,2);
		zC=clump.positions(j,3);
		
		surf(xSph+xC,ySph+yC,zSph+zC,'EdgeColor','none','FaceColor','g','FaceAlpha',1,'FaceColor',color)
	end
end

%% Export clump
% Output is offered in the generic format x_y_z_r. For more specialised
% formats, try the exportClump module.
if isempty(output)==false
	dlmwrite(output, [clump.positions, clump.radii], 'delimiter', ',', 'precision', '%10f')
end

end
