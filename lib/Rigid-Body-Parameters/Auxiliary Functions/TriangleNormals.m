function [FN,TA]=TriangleNormals(TR,flag)
% Calculate normal vectors of individual faces of a triangular surface 
% mesh.
%
% INPUT:
%   - TR    : surface mesh represented as an object of 'TriRep' class,
%             'triangulation' class, or a cell such that TR={Tri,V}, where
%             Tri is an M-by-3 array of faces and V is an N-by-3 array of 
%             vertex coordinates.
%   - flag  : set flag=true to compute normals without using the built-in
%             'faceNormal' function. flag=false is the default setting.
%
% OUTPUT:
%   - FN    : M-by-3 array of triangle normals where, M is the number of 
%             faces. 
%   - TA    : M-by-1 array of triangle areas; obtained as a by-product of 
%             computing FN
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
%


if nargin<2 || isempty(flag), flag=false; end

[Tri,X,fmt]=GetMeshData(TR);

if ~flag && nargout==1
    if fmt>1, TR=triangulation(Tri,X); end
    FN=faceNormal(TR);
    return
end

X1=X(Tri(:,1),:);
X2=X(Tri(:,2),:);
X3=X(Tri(:,3),:);

FN=cross(X2-X1,X3-X1,2);
TA=sqrt(sum(FN.^2,2));
FN=bsxfun(@rdivide,FN,TA);
TA=TA/2;

