function N = patchnormals(FV)
%Vertex normals of a triangulated mesh, area weighted, left-hand-rule
% N = patchnormals(FV) -struct with fields, faces Nx3 and vertices Mx3
%N: vertex normals as Mx3
coord = FV.coord';

%face corners index
A = FV.tri(:,1);
B = FV.tri(:,2);
C = FV.tri(:,3);

%face normals
n = cross(coord(A,:)-coord(B,:),coord(C,:)-coord(A,:)); %area weighted

%vertice normals
N = zeros(size(coord)); %init vertix normals
for i = 1:size(FV.tri,1) %step through faces (a vertex can be reference any number of times)
N(A(i),:) = N(A(i),:)+n(i,:); %sum face normals
N(B(i),:) = N(B(i),:)+n(i,:);
N(C(i),:) = N(C(i),:)+n(i,:);
end
end
