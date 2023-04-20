
function [Z0]= GenerateZ0(mesh, kz, k0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GenerateZ0 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creates Z0 matrix for conditions m>n, m<n, and m=n
%
% inputs
%   mesh        [struct]        generated by onedim.generateMesh
%   kz          [int]      longitudinal wavenumber
%   k0          [int]       wavenumber
% outputs
%   Z0           [N x N double]  impedance matrix
%
% example: 
%
% Roopan -- rroopan@scu.edu
% Logan Barnes-- lbarnes@scu.edu
% Sean Shao -- sshao@scu.edu
% scu 
% 2021
eta=376.73;
A=k0*eta/(2*kz);

% signed difference in index
dM = -(mesh.N-1):1:(mesh.N-1);
% precalculate all exponential terms
EXPTERM = exp(-1j*kz*mesh.dz*abs(dM));
% precalculate one row of the impedance matrix
ZROW = A*(mesh.dz^2)*EXPTERM*onedim.sinc(kz*mesh.dz/(2*pi))^2;
ZROW(mesh.N) = (2*A*mesh.dz/(1j*kz))*(1-exp(-1j*kz*mesh.dz/2)*onedim.sinc(kz*mesh.dz/(2*pi)));

% populate the impedance matrix row by row using a portion of the
% precalculated row.
Z0=zeros(mesh.N,mesh.N);
for m=1:mesh.N
    Z0(m,:) = ZROW((mesh.N-m+1):(end-m+1));
end

end

