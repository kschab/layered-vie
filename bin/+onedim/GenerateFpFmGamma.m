
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GenerateFpFmGamma
% Fp forward going wave operators
% Fm backward going wave operators
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
% inputs
%   mesh           [struct]     generated by onedim.generateMesh(contains epr)
%   hs              [struct]    halfspace structure
%   kz             [float]      longitudinal wavenumber
%   k0              [float]     free space wave number
% outputs
%   Fm             [1 x N double]  Material impedance matrix
%   Fp             [1 x N double]  Material impedance matrix
%
% example: 
%
% Kurt Schab -- kschab@scu.edu
% scu 
% 2021
%Fm- and Fm+ vectors
function [Fm, Fp]= GenerateFpFmGamma(mesh, hs, kz, k0)


% constants
eta=376.73;
kx = sqrt(k0^2 - kz^2);

% halfspace information
eh = hs.er;
h = hs.h;
kh = k0*sqrt(eh);
kzh = sqrt(kh^2-kx^2);
ct = kz/k0;
cth = kzh/(k0*sqrt(eh));
etah = eta/sqrt(eh);

% halfspace reflection coefficient
gh = (etah/cth - eta/ct)/(etah/cth + eta/ct);

% free space contribution
for m=1:mesh.N
    Fm0(m)=-(k0*eta/(2*kz)*mesh.dz)*(exp(-1j*kz*mesh.z(m))*onedim.sinc(kz*mesh.dz/(2*pi)));
    Fp0(m)=-(k0*eta/(2*kz)*mesh.dz)*(exp(1j*kz*mesh.z(m))*onedim.sinc(kz*mesh.dz/(2*pi)));
end

% reflection contribution
for m=1:mesh.N
    FmG(m)=-gh*exp(-1j*kz*2*h)*(k0*eta/(2*kz)*mesh.dz)*(exp(1j*kz*mesh.z(m))*onedim.sinc(kz*mesh.dz/(2*pi)));
end

% construct complete operators
Fp = (1+gh)*Fp0;
Fm = Fm0+FmG;
end