
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GenerateFpFm 
% Fp forward going wave operators
% Fm backward going wave operators
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
% inputs
%   mesh           [struct]     generated by onedim.generateMesh(contains epr)
%   kz             [int]      longitudinal wavenumber
% outputs
%   Fm             [1 x N double]  Material impedance matrix
%   Fp             [1 x N double]  Material impedance matrix
%
% example: 
%
% Roopan -- rroopan@scu.edu
% Logan Barnes-- lbarnes@scu.edu
% Sean Shao -- sshao@scu.edu
% scu 
% 2021
%Fm- and Fm+ vectors
function [Fm, Fp]= GenerateFpFm(mesh, kz, k0)
eta=376.73;
for m=1:mesh.N
    Fm(m)=-(k0*eta/(2*kz)*mesh.dz)*(exp(-1j*kz*mesh.z(m))*onedim.sinc(kz*mesh.dz/(2*pi)));
    Fp(m)=-(k0*eta/(2*kz)*mesh.dz)*(exp(1j*kz*mesh.z(m))*onedim.sinc(kz*mesh.dz/(2*pi)));
end
end