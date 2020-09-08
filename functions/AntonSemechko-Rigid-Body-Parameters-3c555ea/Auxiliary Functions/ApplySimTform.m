function X=ApplySimTform(X,T,opt)
% Apply linear transformation to a mesh or a set of 3D points.
%
%   - X     : mesh or N-by-3 list of Cartesian point coordinates
%   - T     : 3-by-3 rotation matrix or 4-by-4 homogeneous transformation
%             matrix
%   - opt   : if opt=true then output equals T*X {default}, and 
%             T\X otherwise.
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
%


if nargin<2 || isempty(X) || isempty(T)
    error('Insufficient number of input arguments')
end

d=size(T,1);
if ~isnumeric(T) || ~ismatrix(T) || d~=size(T,2) || d<3 || d>4 || (d==4 && ~isequal(T(4,:),[0 0 0 1]))
    error('Invalid format for 2nd input argument')
end
    
if nargin<3 || isempty(opt), opt=true; end

% Face and vertex lists
fmt=0;
if ~(isnumeric(X) && ismatrix(X) && size(X,2)==3)
    [Tri,X,fmt]=GetMeshData(X);
end

% Apply transformation
if d==4, X(:,4)=1; end
if opt % "forward" transformation
    X=(T*(X'))';
else % inverse transformation
    X=(T\(X'))';
end
X=X(:,1:3);

% Output 
switch fmt
    case 1
        X=triangulation(Tri,X);
    case 2
        X=TriRep(Tri,X); %#ok<*DTRIREP>
    case 3
        X={Tri X};
    case 4
        X=struct('faces',Tri,'vertices',X);
end

