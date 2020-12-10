function [h,RBP]=VisualizeLocalFrame(TR,prm,ax_dir)
% Visualize local frame of reference of an objected represented by a
% closed, triangular surface mesh.
%
% INPUT:
%   - TR     : triangular surface mesh specified in one of the following
%              formats:
%               a. 'TriRep' object (may be obsolete depending on your 
%                   current version of Matlab)
%               b. 'triangulation' object
%               c. Structure with the exact same fields as the one returned 
%                   by the 'isosurface' function; that is  
%                       TR.faces    : M-by-3 list of faces
%                       TR.vertices : N-by-3 list of vertex co-ordinates
%               d. Cell such TR={Tri X} where Tri is a M-by-3 list of faces
%                  and X is a N-by-3 list of vertex co-ordinates
%   - prm    : integer specifying type of primitive to use for coordinate
%              frame visualization. The three available options are:
%               1. Bounding-box <-- {default}
%               2. Approximating cuboid
%               3. Approximating ellipsoid
%   - ax_dir : optional 1-by-3 vector, [sx sy sz], used to determine 
%              positive senses of the prinicpal axes. For example, 
%              setting ax_dir=[-1 -1 1] will cause +ve directions of the
%              1st and 2nd primary axes to flip. ax_dir=[1 1 1] is the 
%              default setting. No that to ensure that visualized 
%              coordinate frame maintains right-hand convention verify that
%              prod(ax_dir)==1.
%
% OUTPUT:
%   - h     : vector containing handles for all graphical objects, so that
%             h=[h_prm h_arw] where h_prm are the handles returned by the 
%             the 'RBP_cuboid' or 'RBP_ellipsoid' functions and h_arw
%             contains handles for the arrows used to represent +ve
%             directions of the principal axes.            .
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
%


if nargin<2 || isempty(prm)
    prm=1;
elseif numel(prm)~=1 || ~isnumeric(prm) || sum((1:3)==prm)==0
    error('Invalid entry for 2nd input argument (''prm'')')
end

if nargin<3 || isempty(ax_dir)
    ax_dir=ones(1,3);
elseif numel(ax_dir)~=3 || ~isvector(ax_dir) || sum(ax_dir==1 |ax_dir==-1)~=3
    error('Invalid entry for 3rd input argument (''ax_dir'')')
end

% Inertia-based reference frame
[R,~,~,RBP,BB,bb_prm]=InertiaBasedLocalFrame(TR);
R=bsxfun(@times,R,ax_dir(:)');

% Visualize bounding-box or one of approximating primitives
switch prm
    case 1 % bounding-box
        
        abc=bb_prm.abc;        
        [Tri,V]=GetMeshData(TR);
        
        figure('color','w')
        h1=patch('faces',Tri,'vertices',V);
        set(h1,'EdgeColor','none','FaceColor',0.4*ones(3,1),'EdgeLighting','none')
        hold on
        
        h2=patch(BB);
        col=[0.85,0.33,0.10];
        set(h2,'EdgeColor','k','FaceColor',col,'FaceAlpha',0.1,'LineWidth',2)
        
        hL1=camlight('headlight');
        set(hL1,'style','infinite','position',bb_prm.centroid+20*max(abc)*R(:,1)')
        hL2=light('position',bb_prm.centroid-20*max(abc)*R(:,1)');
        set(hL2,'style','infinite')
        lighting phong
        
        axis equal off %vis3d
        
        h_prm=[h1 h2 hL1 hL2];

        abc=abc/2;
        
    case 2 % approximating cuboid
        
        [RBP,abc,~,h_prm]=RBP_cuboid(TR,true);
        abc=abc/2;
        
    case 3 % approxiamting ellipsoid
        
        [RBP,abc,~,h_prm]=RBP_ellipsoid(TR,true);
        
end
set(h_prm(1),'FaceAlpha',0.5,'FaceColor',0.5*[1 1 1],'EdgeColor','none')

pL1=get(h_prm(end-1),'position');
pL2=get(h_prm(end)  ,'position');
delete(h_prm((end-1):end)) 

% Visualize local frame 
if prm==1
    C=bb_prm.centroid;
else
    C=RBP.centroid;
end
abc(:)=min(abc);
col={'r' 'g' 'b'};
r=abc(3)/15;          % arrow radius
hf=1-0.4*abc(3)./abc; % arrow head fraction
h_arw=[0 0 0];
for j=1:3
    Vj=R(:,j);
    hj=arrow3d(abc(j)*[0 Vj(1)]+C(1),abc(j)*[0 Vj(2)]+C(2),abc(j)*[0 Vj(3)]+C(3),hf(j),r,2*r,col{j});
    h_arw(j)=hj;
end

hL1=light; set(hL1,'style','infinite','position',pL1)
hL2=light; set(hL2,'style','infinite','position',pL2)
lighting phong
h_prm(end-1)=hL1; 
h_prm(end)  =hL2; 

if nargout>0, h=[h_prm h_arw]; end

