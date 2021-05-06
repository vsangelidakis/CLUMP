function [RBP,abc,R,h]=RBP_ellipsoid(TR,vis,varargin)
% Compute rigid-body parameters of an object represented by a triangular 
% surface mesh and then find an approximating ellipsoid that has the same
% inertial parameters. Note that while the approximating ellipsoid will 
% have the same inertial parameters as the input object, in general, its 
% volume will NOT be the same as that of the input object.
%
% INPUT:
%   - TR    : triangular surface mesh specified in one of the following
%             formats:
%               a. 'TriRep' object (may be obsolete depending on your 
%                   current version of Matlab)
%               b. 'triangulation' object
%               c. Structure with the exact same fields as the one returned 
%                   by the 'isosurface' function; that is  
%                       TR.faces    : M-by-3 list of faces
%                       TR.vertices : N-by-3 list of vertex co-ordinates
%               d. Cell such TR={Tri X} where Tri is a M-by-3 list of faces
%                  and X is a N-by-3 list of vertex co-ordinates
%   - vis   : optional binary input argument. vis=true if you want to
%             visualize the input mesh along with the approximating 
%             ellipsoid and vis=false otherwise. The latter is the default 
%             setting.
%
% OUTPUT:
%   - RBP   : structure containin rigid-body parameters of the object 
%             represented by TR. See 'RigidBodyParams' function for more 
%             info. 
%   - abc   : 1-by-3 vector containing half-lengths of the approximating 
%             ellipsoid's primary axes in order of decreasing magnitude.
%   - R     : 3-by-3 rotation matrix such that R(i,:) is the direction of 
%             the i-th primary axis. Note that R(i,:) corresponds to abc(i). 
%   - h     : vector containing graphical object handles, 
%             h=[h1 h2 ha hb hc hL1 hL2], where: 
%               - h1        : handle for the patch representing input mesh
%               - h2        : handle for patch representing the primitive
%               - ha,hb,hc  : line handles of primary ellipses
%               - hL1,hL2   : lighting objects 
%           
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
%


h=[];
if nargin<2 || isempty(vis), vis=false; end
vis=vis(1);
if ~vis && nargout<1, return; end

% Find RBPs 
if nargin<3 || isempty(varargin{1})
    RBP=RigidBodyParams(TR);
else
    RBP=varargin{1};
end

% Find ellipsoid parameters
L=RBP.eigs;
l1=L(1)-L(3)+L(2);
l2=L(1)-L(2)+L(3);
l3=L(3)-L(1)+L(2);

g=1/(8*pi/15);
a=(g*l1^2/sqrt(l2*l3))^(1/5);
b=(g*l2^2/sqrt(l1*l3))^(1/5);
c=(g*l3^2/sqrt(l1*l2))^(1/5);
abc=[a b c];

% Local frame of reference
R=fliplr(RBP.PAI);
R(:,3)=cross(R(:,1),R(:,2));

% Construct a mesh for the ellipsoid. Note that this mesh is only an
% approximation of the true ellipsoid, and therefore its inertial 
% parameters will approach those of TR in the limit N-->Inf, where N is 
% the number of mesh vertices.
if ~vis(1), return; end
tr=SubdivideSphericalMesh(IcosahedronMesh,5);
X=tr.Points;
X(:,1)=a*X(:,1);
X(:,2)=b*X(:,2);
X(:,3)=c*X(:,3);
X=bsxfun(@plus,(R*X')',RBP.centroid);

Tri=tr.ConnectivityList;
tr=struct();
tr.faces=Tri;
tr.vertices=X;

% Construct great ellipses
t=linspace(0,2*pi,1E3);
x=cos(t);
y=sin(t);

[GC1,GC2,GC3]=deal(zeros(1E3,3));

GC1(:,1)=a*x; GC1(:,2)=b*y; % xy-plane
GC2(:,1)=a*y; GC2(:,3)=c*x; % zx-plane
GC3(:,2)=b*x; GC3(:,3)=c*y; % yz-plane

C=RBP.centroid;
GC1=bsxfun(@plus,(R*GC1')',C);
GC2=bsxfun(@plus,(R*GC2')',C);
GC3=bsxfun(@plus,(R*GC3')',C);

% Visualize meshes
[Tri,V]=GetMeshData(TR);

figure('color','w')
h1=patch('faces',Tri,'vertices',V); 
set(h1,'EdgeColor','none','FaceColor',0.4*ones(3,1),'EdgeLighting','none')
hold on

h2=patch(tr); 
col=[0.85,0.33,0.10];
set(h2,'EdgeColor','none','FaceColor',col,'FaceAlpha',0.1)

hL1=camlight('headlight');
set(hL1,'style','infinite','position',RBP.centroid+20*max(abc)*R(:,1)')
hL2=light('position',RBP.centroid-20*max(abc)*R(:,1)');
set(hL2,'style','infinite')
lighting phong

ha=plot3(GC1(:,1),GC1(:,2),GC1(:,3),'-k','LineWidth',2);
hb=plot3(GC2(:,1),GC2(:,2),GC2(:,3),'-k','LineWidth',2);
hc=plot3(GC3(:,1),GC3(:,2),GC3(:,3),'-k','LineWidth',2);
axis equal off vis3d

h=[h1 h2 ha hb hc hL1 hL2];
