function fv=unit_cube_mesh
% Construct quadrilateral mesh of a unit cube with edges along
% positive x-, y-, and z-axes, and one corner at the origin. 
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
%

X=[0 0 0; ...
   1 0 0; ...
   1 1 0; ...
   0 1 0; ...
   0 0 1; ...
   1 0 1; ...
   1 1 1; ...
   0 1 1];

F=[1 4 3 2;
   5 6 7 8;
   2 3 7 6;
   3 4 8 7;
   1 5 8 4;
   1 2 6 5];
 
fv.faces=F;
fv.vertices=X;
 
