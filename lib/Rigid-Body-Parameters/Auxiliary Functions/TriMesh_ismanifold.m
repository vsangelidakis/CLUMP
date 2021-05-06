function [flag,E,FB,EA,Ne]=TriMesh_ismanifold(TR)
% Check if triangular surface mesh is manifold; that is every non-boundary 
% edge is shared by two faces.
%
% INPUT:
%   - TR  : triangular surface mesh represented as an object of 'TriRep' 
%           class, 'triangulation' class, or a cell such that TR={Tri,V},
%           where Tri is an M-by-3 array of faces and V is an N-by-3 array
%           of vertex coordinates.
%
% OUTPUT:
%   - flag : true if TR is manifold and false otherwise.
%   - E    : non-boundary edge array
%   - FB   : boundary edge array
%   - EA   : non-boundary edge attachment cell 
%   - Ne   : number of elements in each EA cell
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
%


[Tri,V,fmt]=GetMeshData(TR);
if size(Tri,2)~=3
    error('This function is intended for TRIANGULAR surface meshes.')
end

%warning('off','MATLAB:triangulation:PtsNotInTriWarnId')
if fmt>1
    TR=triangulation(Tri,V);
end

% Does TR have a boundary?
E=sort(edges(TR),2);
FB=freeBoundary(TR);

flag_closed=isempty(FB);
if ~flag_closed
    FB=sort(FB,2);
    idx=ismember(E,FB,'rows');
    E(idx,:)=[];   
end

% Edge attachments
EA=edgeAttachments(TR,E);
Ne=cellfun(@length,EA);
flag=all(Ne==2);
if ~flag && ~flag_closed
    EAb=edgeAttachments(TR,FB);
    Neb=cellfun(@length,EAb);
    flag=all(Neb==1);
end
%warning('on','MATLAB:triangulation:PtsNotInTriWarnId')
