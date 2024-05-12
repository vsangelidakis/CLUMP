function [TR,cnt]=ConsistentNormalOrientation(TR,ord)
% Make sure all faces of a triangular surface mesh have consistent and
% proper orientation. Proper orientation is one where face normals point 
% outside the region enclosed by a mesh. 
%
% INPUT:
%   - TR  : closed or open triangular surface mesh represented as an object
%           of 'TriRep' class, 'triangulation' class, or a cell such that
%           TR={Tri,V}, where Tri is an M-by-3 array of faces and V is an 
%           N-by-3 array of vertex coordinates.
%   - ord : (optional) set opt=true to enforce orignal order of faces in 
%           the face-vertex connectivty list. If the relative order of 
%           faces is not important, use the dafault opt=false setting to 
%           save time. 
%
% OUTPUT:
%   - TR  : mesh with consistently oriented faces. Format is the same as 
%           that of the input. Note, relative position of faces in the 
%           face-vertex connectivity list of the output mesh may differ 
%           from the input even though none of the faces were inverted. 
%           To enforce original order make sure to set the option ord=true.
%   - cnt : number of inverted faces.
%   
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
%


if nargin<1 || isempty(TR)
    error('Insufficient number of input arguments')
end

if nargin<2 || isempty(ord)
    ord=nargout>1;
elseif ~islogical(ord) || numel(ord)~=1
    error('Invalid entry for 2nd input argument (ord)')
end


% Get mesh data
[Tri,V,fmt]=GetMeshData(TR);
if size(Tri,2)~=3
    error('This function is intended for TRIANGULAR surface meshes.')
end

if fmt>1
    TR=triangulation(Tri,V);
end


% Is the mesh manifold?
cnt=0;
[chk_mnfd,E,FB,EAc]=TriMesh_ismanifold(TR);
if ~chk_mnfd
    error('Mesh is non-manifold')
end
chk_clsd=isempty(FB);


% Quick test for potentially inverted faces using dihedral angles
for i=1:size(E,1), E(i,:)=EAc{i}; end
FN=faceNormal(TR);
DA=sum(FN(E(:,1),:).*FN(E(:,2),:),2);
chk=DA<0;
flag=any(chk);


% "Feasible" edge
FE=[Tri(1,:);circshift(Tri(1,:),[0 -1])]';  

Nf=size(Tri,1);
TriNew=zeros(Nf,3);
TriNew(1,:)=Tri(1,:);
idx_pass=false(Nf,1);
idx_pass(1)=true;


% Connectivity list
% -------------------------------------------------------------------------
i=1; 
%pc=round((0.05:0.05:1)*Nf);
%fprintf('%-4u',5:5:100)
%fprintf('%%\n')

f=zeros(1,3);
idx2=false(1,3);
idx_c=cell(1);
while i<Nf && flag    
    
    % Get face attached to the first edge in the FE list    
    if isempty(FE)
        FE=get_feasible_edges(TriNew(1:i,:),V);
        if isempty(FE)
            error('Face-vertex connectivty list either contains duplicate faces or is a union of two or more disconnected triangle sets.')
        end
    end

    idx_c(:)=edgeAttachments(TR,FE(1,:));
    fa=idx_c{1};
    fa(idx_pass(fa))=[]; % remove faces that were already processed
    if isempty(fa)
        FE(1,:)=[];
        continue
    end
    fa=fa(1);
    idx_pass(fa)=true;
    
    i=i+1; % found unprocessed face
    
    % New vertex sequence
    f(:)=Tri(fa,:);
    idx2(:)=f==FE(1,1) | f==FE(1,2);
    f(:)=[FE(1,2) FE(1,1) f(~idx2)]; 
    TriNew(i,:)=f;
    
    % Add more feasible edges
    FEi=[FE(1,1) f(3); f(3) FE(1,2)];
    FE=[FE;FEi];     %#ok<*AGROW>
    FE(1,:)=[];
    
    %if i==pc(1)
    %    fprintf('.   ');
    %    pc(1)=[];
    %end
    
end
   

% Decide on global orientation of the normals
% -------------------------------------------------------------------------
if ~flag, TriNew=Tri; end

if chk_clsd % mesh is closed
    
    Vol=ClosedMeshVolume({TriNew V});
    if Vol<0
        TriNew=fliplr(TriNew);
    end
    
else % mesh is open, need a different approach to determine normal orientations
    
    [FN,W]=TriangleNormals({TriNew V},true);
    
    W=W/sum(W);
    Vf=(V(TriNew(:,1),:)+V(TriNew(:,2),:)+V(TriNew(:,3),:))/3;
    C=W'*Vf;
    
    D=bsxfun(@minus,Vf,C);
    %D=bsxfun(@rdivide,D,sqrt(sum(D.^2,2))+eps);
    
    s=sum(FN.*D,2)>0;
    if sum(W(s))<0.5, TriNew=fliplr(TriNew); end
    
end


% Enforce original face order
% -------------------------------------------------------------------------
if ord
    
    % Map new faces to the old ones
    if flag    
        [~,id_fwd]=ismember(sort(TriNew,2),sort(Tri,2),'rows');
        [~,id_srt]=sort(id_fwd);
        TriNew=TriNew(id_srt,:);
    end
        
    % Match first vertex of each face
    for i=1:2
        idx=Tri(:,1)==TriNew(:,i+1);
        if any(idx)
            tri=TriNew(idx,:);
            tri=circshift(tri,[0 -i]);
            TriNew(idx,:)=tri;
        end
    end
           
    % Count the number of inverted faces
    if nargout>1
        cnt=nnz(sum(Tri==TriNew,2)<3);
    end
    
end
clear Tri


switch fmt
    case 1
        TR=triangulation(TriNew,V);
    case 2
        TR=TriRep(TriNew,V); %#ok<*DTRIREP>
    case 3
        TR={TriNew V};
    case 4
        TR=struct('faces',TriNew,'vertices',V);
end



function FE=get_feasible_edges(Tri,V)
% Find boundary edges of triangle set that has been processed.

% Remove non-referenced vertices
[idx_v,~,idx]=unique(Tri(:));

Nv_unq=length(idx_v);
v_id=(1:Nv_unq)';
Tri2=v_id(idx);
Tri2=reshape(Tri2,[],3);
TR2=triangulation(Tri2,V(idx_v,:));

% Boundary edges
FE=[];
BE=freeBoundary(TR2);
if isempty(BE), return; end

% BEs as they appear in Tri
chk=ismember(Tri2,BE(:));

E=cat(1,Tri(:,[1 2]), Tri(:,[2 3]), Tri(:,[3 1]));
idx=cat(1,chk(:,1) & chk(:,2),chk(:,2) & chk(:,3),chk(:,3) & chk(:,1));

FE=E(idx,:);

