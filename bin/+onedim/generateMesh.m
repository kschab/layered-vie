function [mesh]= generateMesh (N, h)
%       output      name        input
mesh.N=N;
mesh.h=h;
mesh.dz=h/N;
mesh.z=([1:mesh.N]-.5)*mesh.dz;

end