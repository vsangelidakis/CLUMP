function [RBP,abc,R,h]=RBP_cuboid(TR,vis,varargin)
% Compute rigid-body parameters of an object represented by a triangular 
% surface mesh and then find an approximating cuboid that has the same
% inertial parameters. Note that while the approximating cuboid will 
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
%             cuboid and vis=false otherwise. The latter is the default 
%             setting.
%
% OUTPUT:
%   - RBP   : structure containing rigid-body parameters of the object 
%             represented by TR. See 'RigidBodyParams' function for more 
%             info. 
%   - abc   : 1-by-3 vector containing length, width and height parameters 
%			  of the approximating cuboid.  
%   - R     : 3-by-3 rotation matrix such that R(i,:) is the direction of 
%             the i-th dimension, so that R(i,:) corresponds to abc(i). 
%   - h     : vector containing graphical object handles, h=[h1 h2 hL1 hL2], 
%             where: 
%               - h1       : handle for the patch representing input mesh
%               - h2       : handle for patch representing the primitive
%               - hL1,hL2  : lighting objects       
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

% Find cuboid parameters
L=RBP.eigs;
l1=L(1)-L(3)+L(2);
l2=L(1)-L(2)+L(3);
l3=L(3)-L(1)+L(2);

a=(6*l1^2/sqrt(l2*l3))^(1/5);
b=(6*l2^2/sqrt(l1*l3))^(1/5);
c=(6*l3^2/sqrt(l1*l2))^(1/5);
abc=[a b c];

% Local frame of reference
R=fliplr(RBP.PAI);
R(:,3)=cross(R(:,1),R(:,2));

% Construct a quad mesh for the cuboid
if ~vis(1), return; end
X=[1 1; -1 1; -1 -1; 1 -1];
X=[X;X];
X(:,3)=1;
X(5:8,3)=-1;

X(:,1)=(a/2)*X(:,1);
X(:,2)=(b/2)*X(:,2);
X(:,3)=(c/2)*X(:,3);
X=bsxfun(@plus,(R*X')',RBP.centroid);

F=[1 2 3 4;...
   5 8 7 6;...
   8 5 1 4;
   2 6 7 3;
   1 5 6 2;
   3 7 8 4];

tr.faces=F;
tr.vertices=X;

% Visualize meshes
[Tri,V]=GetMeshData(TR);

figure('color','w')
h1=patch('faces',Tri,'vertices',V); 
set(h1,'EdgeColor','none','FaceColor',0.4*ones(3,1),'EdgeLighting','none')
hold on

h2=patch(tr); 
col=[0.85,0.33,0.10];
set(h2,'EdgeColor','k','FaceColor',col,'FaceAlpha',0.1,'LineWidth',2)

hL1=camlight('headlight');
set(hL1,'style','infinite','position',RBP.centroid+40*max(abc)*R(:,1)')
hL2=light('position',RBP.centroid-40*max(abc)*R(:,1)');
set(hL2,'style','infinite')
lighting phong

axis equal off vis3d

h=[h1 h2 hL1 hL2];

