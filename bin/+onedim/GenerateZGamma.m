
function [ZG]= GenerateZGamma(mesh, hs, kz, k0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GenerateZGamma -- HALFSPACE VERSION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creates reflection impedance matrix Z^\Gamma -- HALFSPACE VERSION
%
% inputs
%   mesh        [struct]        generated by onedim.generateMesh
%   hs          [struct]   halfspace object
%   kz          [int]      longitudinal wavenumber
%   k0          [int]       wavenumber
% outputs
%   Z           [N x N double]  reflection impedance matrix
%
% example: 
%
% Kurt Schab -- kschab@scu.edu
% scu 
% 2021

% constants
eta=376.73;
A=k0*eta/(2*kz);
kx = sqrt(k0^2 - kz^2);
dz = mesh.dz;

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

% construct reflection impedance matrix
ZZ = gh*dz^2*A*onedim.sinc(kz*dz/2)*exp(-1j*kz*2*h);
ZG = zeros(mesh.N,mesh.N);
for m = 1:mesh.N
    EM = exp(1j*kz*mesh.z(m));
    for n = 1:mesh.N
        EN = exp(1j*kz*mesh.z(n));
        ZG(m,n) = ZZ*EM*EN;
    end
end
        

end
