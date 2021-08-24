function [faces,vertices]=ExtractSurface(clump, N_sphere, N_circle, visualise)
%% Tesselation of the surface of a clump into a surface mesh
% 2021 © V. Angelidakis, S. Nadimi, M. Otsubo, S. Utili.

%% INPUT:
%	- clump    : either "clump" object or N x 4 matrix with columns of [x,y,z,r], where x,y,z the centroid of each sphere and r its radius
%	- N_sphere : Number of vertices on the surface of each member-sphere of the clump
%	- N_circle : Number of vertices on the circle defined as the intersection of two overlapping spheres
%	- visualise: Boolean whether to plot the generated surface mesh of the clump surface

%% OUTPUT:
%	- faces	   : faces of generated surface mesh
%	- vertices : vertices of generated surface mesh

%% EXAMPLE
%
% N_sphere=400;
% N_circle=200;
% clump=[
% 	1,0,0,1.1;
% 	2,1,0,1.1;
% 	3,0,0,1.2;
%	];
% [faces,vertices]=ExtractSurface(clump,N_sphere,N_circle,visualise)

%% Check format of input
switch class(clump)
	case 'struct'
		if isfield(clump,'positions') && isfield(clump,'radii')
			if size(clump.positions,2)<3
				error('Invalid format; clump.positions should have size N x 3!')
			end
			if size(clump.radii,2)>1
				error('Invalid format; clump.radii should have size N x 1!')
			end
			spheresList=[clump.positions,clump.radii]; % TODO: Remove the duplicate variable later on
			
		else
			error('Invalid format! The struct should have fields "positions" and "radii"!')
		end
	case 'double'
		if size(clump,2)~=4
			error('Invalid format, should be x,y,z,r!')
		end
		spheresList=clump; % TODO: Remove the duplicate variable later on
end

[x,y,z,r]=deal(spheresList(:,1),spheresList(:,2),spheresList(:,3),spheresList(:,4));

%% Main body of the function
%% Import Dependencies
addpath(genpath('MyCrust'))

%% Contact detection between all spheres (all possible combinations) - Record interactions
interactions=[];
ind=1;
for i=1:size(spheresList,1)-1
	for j=i+1:size(spheresList,1)
		if i==j
			continue
		end
		inContact=sphereContact(spheresList(i,:),spheresList(j,:));
		if inContact
			interactions(ind,1:2)=[i,j];
			ind=ind+1;
		end
	end
end

%% Generate points for each sphere
for i=1:size(spheresList,1)
	[S{i}.vertices,S{i}.faces]=makeSphere(x(i),y(i),z(i),r(i),N_sphere);
end

%% Perform contact detection to detect and delete points of each sphere that are included in other spheres, in order to get only the points on the surface of the clump
for i=1:size(interactions,1)
	
	% For interaction [sphere1,sphere2], check which vertices of sphere1 are inside sphere2
	for j=size(S{1,interactions(i,1)}.vertices,1):-1:1 % start deleting from end to start
		if spherePotential(S{1,interactions(i,1)}.vertices(j,:),spheresList(interactions(i,2),:),true)
			v=S{1,interactions(i,1)}.vertices(j,:);
			% 			scatter3(v(:,1),v(:,2),v(:,3),20,'b','filled')
			S{1,interactions(i,1)}.vertices(j,:)=[];
		end
	end
	
	% For interaction [sphere1,sphere2], check which vertices of sphere2 are inside sphere1
	for j=size(S{1,interactions(i,2)}.vertices,1):-1:1 % start deleting from end to start
		if spherePotential(S{1,interactions(i,2)}.vertices(j,:),spheresList(interactions(i,1),:),true)
			v=S{1,interactions(i,2)}.vertices(j,:);
			% 			scatter3(v(:,1),v(:,2),v(:,3),20,'b','filled')
			S{1,interactions(i,2)}.vertices(j,:)=[];
		end
	end
	
end
% scatter3(S{1,1}.vertices(:,1),S{1,1}.vertices(:,2),S{1,1}.vertices(:,3),'filled')
% scatter3(S{1,2}.vertices(:,1),S{1,2}.vertices(:,2),S{1,2}.vertices(:,3),'filled')


%% Calculate points on the intersection of each pair of interacting spheres
vertices=[];
for i=1:size(interactions,1)
	n=spheresList(interactions(i,2),1:3)-spheresList(interactions(i,1),1:3); % (not normalised) normal vector of each interaction
	
	d=norm(n); % centroidal distance between sphere1-sphere2 in each interaction
	n=n/norm(n); % normalised normal vector of each interaction
	
	r1=spheresList(interactions(i,1),4); % radius of sphere1
	r2=spheresList(interactions(i,2),4); % radius of sphere2
	
	h = sqrt( (2*r1*d)^2 - (r1^2 + d^2 - r2^2)^2 )/(2*d); % Radius of intersection circle
	alph=acos( (r1^1+d^2-r2^2) / (2*r1*d) );
	h1=r1*(1-cos(alph));
	
	C=spheresList(interactions(i,1),1:3)+n*(r1-h1); % Contact point
	
	% FIXME: Revisit the calculation of the perpendicular vector
	n3=n;
	n1=[n(3) 0 -n(1)]; % Vector perpendicular to n
	if norm(n1)==0
		n1=[n(2) 0 -n(1)];
	end
	n1=n1/norm(n1); % Normalise n1
	n2=cross(n3,n1);
	% 	dot(n1,n3)
	
	% Generate points of intersection circle
	a=-2*pi:pi/(N_circle/4):2*pi;
	% For each circle point.
	px = C(1) + h * (n1(1) * cos(a) + n2(1) * sin(a));
	py = C(2) + h * (n1(2) * cos(a) + n2(2) * sin(a));
	pz = C(3) + h * (n1(3) * cos(a) + n2(3) * sin(a));
	
	if imag(px(1,1))>0
		% 		i
		break
	end
	
	%% TODO: Turn px, py, pz calculation into matrix operation
	
	% 	if ~isnan(px(1,1))
	S{1,interactions(i,1)}.circlevertices=[px' py' pz']; %% FIXME: This line is wrong! S{1,i} instead?
	% 	S{1,i}.circlevertices=[px' py' pz']; %% FIXME: This line is wrong! S{1,i} instead?
	
	vertices=[vertices;[px' py' pz']];
	% 		scatter3(px,py,pz,20,'filled','r')
	% 	end
	
	%% TODO: Instead of adding all vertices apriori, I could first check that the new vertices are not inside existing spheres. They can be on (i.e. potential=0) but not zero. To this, I need a new contact detection, allowing zero
	
end

%% Collect vertices from all spheres in one variable
for i=1:size(S,2)
	vertices=[vertices;S{1,i}.vertices];
end
vertices = real(unique(vertices,'rows')); % FIXME: Come back to check when we need the real() function

%% Generate mesh using the Crust algorithm (Amenta et al, 1999)
% p=vertices; %% TODO: Remove the duplication: Keep only vertices variable
[faces,~]=MyRobustCrust(vertices);
faces=double(faces); % transform from int32 to double

%% TODO: Improve the case where some of the spheres are not in contact with the rest of the clump

if visualise
	fig=figure('Position',[200 200 600 600]);
	box on; grid on; hold on;
	axis vis3d equal
	patch('Faces',faces,'vertices',vertices,'FaceColor','c','EdgeColor','none')
	% 	trisurf(faces,vertices(:,1),vertices(:,2),vertices(:,3),'facecolor','c','edgecolor','none')
	alpha 0.5
	view(3)
	set(gca,'visible','off')
	camlight
	camproj('perspective')
end
end


function inContact=sphereContact(sphere1,sphere2)
%% Function to perform contact detection between two spheres
% inContact: boolean: whether sphere1 and sphere2 intersect
% sphere1:	[1 x 4] [x,y,z,r]:	test sphere 1
% sphere2:	[1 x 4] [x,y,z,r]:	test sphere 2

d0=norm(sphere2(1:3)-sphere1(1:3)); % Centroidal distance of the spheres
if d0<=(sphere1(4)+sphere2(4))
	inContact=true;
else
	inContact=false;
end
% inContact=sqrt( ( (sphere(1)-point(1))^2 + (sphere(2)-point(2))^2 + (sphere(3)-point(3))^2 )/(sphere(4))^2 ) - 1 <= 0;
end


function isInside=spherePotential(point,sphere,allowZero)
%% Function to determine whether a point is inside a sphere
% isInside: boolean: whether the test point is inside the sphere of interest
% point:	[1 x 3] x,y,z:	test point
% sphere:	[1 x 4] x,y,z,r	sphere of interest
% allowZero: boolean: whether to consider 0 values as contact, i.e. returning true
if allowZero
	isInside=sqrt( ( (sphere(1)-point(1))^2 + (sphere(2)-point(2))^2 + (sphere(3)-point(3))^2 )/(sphere(4))^2 ) - 1 <= 0;
else
	isInside=sqrt( ( (sphere(1)-point(1))^2 + (sphere(2)-point(2))^2 + (sphere(3)-point(3))^2 )/(sphere(4))^2 ) - 1 < 0;
end
end


function [vertices,faces]=makeSphere(X,Y,Z,radius,N) %radius
%% Function to create a surface mesh of a sphere with radius r, centered at (x,y,z) with N vertices.
% Returns vertices/faces

vertices=zeros(N,3);
inc=pi*(3-sqrt(5));
off=2/N;
for k=0:N-1
	y=k*off-1+off/2;
	r=sqrt(1-y^2);
	phi=k*inc;
	vertices(k+1,1:3)=[cos(phi)*r*radius, y*radius, sin(phi)*r*radius];
end
vertices=vertices+[X,Y,Z];

faces=convhull(vertices);
% patch('vertices',vertices,'Faces',k,'FaceColor','g')
end