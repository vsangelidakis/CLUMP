function TA=TriangleAreas(TR)
% Calculate area of individual faces of a triangular surface mesh.
%
% INPUT:
%   - TR    : surface  mesh represented as an object of 'TriRep' class,
%             'triangulation' class, or a cell such that TR={Tri,V}, where
%             Tri is an M-by-3 array of faces and V is an N-by-3 array of 
%             vertex coordinates.
%
% OUTPUT:
%   - TA    : M-by-1 array of triangle areas where, M is the number of 
%             faces. 
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
%


% Face and vertex lists
[Tri,V]=GetMeshData(TR);
if size(Tri,2)~=3
    error('This function is meant ONLY for triangular surface meshes')
end

% Coordinates of triangles' vertices
V1=V(Tri(:,1),:);
V2=V(Tri(:,2),:);
V3=V(Tri(:,3),:);

% Direction vectors
D1=V2-V1;
D2=V3-V1;

% Triangle areas
TA=cross(D1,D2,2);
TA=sum(TA.^2,2);
TA=sqrt(TA)/2;

