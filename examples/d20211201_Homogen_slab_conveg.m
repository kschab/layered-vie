close all
clear
clc

% template for a ``complete'' scattering analysis using the onedim library
% for one-dimensional method of moments problems.
%
% think of this program as a do-it-all demo of the library's capabilities. 
% practical analyses will use bits and pieces of this code as needed.
%
% scu em lab
% 2021
 e0_ = [3-1j*0.2 6-1j*0.2 5 1j*0.20 -1j*100]
for e0 = e0_ 
%% I: user defined parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_ = floor(logspace(1,4,25));
Error_ = N_*0;
index = 0;
for N = N_
   N
    index = index+1;
% excitation parameters
thetai = pi/3;
k0 =  1;
mode = 'homogeneous';
% e0 = 3-1j*0.2;     % just a parameter defining arbitrary permittivity profiles

% geometry parameters
k0h = 6.4;  % free space electrical thickness of slab

k0zext = [-k0h, 2*k0h]; % free space electrical size of extended mesh for plotting 
Next = 75; % number of points in each exterior plot region

% permittivity profile -- demo of two ways of defining the layered media
switch mode
    case 'function'     % define the permittivity based on a function
        epr = @(z) z*e0;
    case 'mask'         % define the permittivity based on a binary mask
        mask = randi(2,N,1)-1;
        epr = ones(N,1)*1;
        epr(mask>0) = e0;
    case 'homogeneous'  % homogeneous slab
        epr = ones(N,1)*e0;
end

%% II: setup 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create mesh
h = k0h/k0;
mesh = onedim.generateMesh(N,h);

% assign permittivity (two options given here)
switch mode
    case 'function'
        mesh.epr = onedim.assignPermFun(mesh,epr);
    otherwise
        mesh.epr = epr;
end

% compute transverse and longitudinal wavenumebers
kx = k0*sin(thetai);
k0z = k0*sqrt(1-sin(thetai)^2);

%% III: calculate impedance matrices, excitation vector, + other operators
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute incident field
V = onedim.incidentPlusPlaneWave(mesh,k0z);

% compute vacuum impedance matrix 
Z0 = onedim.GenerateZ0(mesh, k0z, k0);

% compute material matrix 
Zp = onedim.GenerateZp(mesh,k0); %convert resistivity

% compute forward and backward going wave operators 
[Fm, Fp]= onedim.GenerateFpFm(mesh, k0z, k0);

%% IV: calculate driven solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute currents
I = (Z0+Zp)\V;

% compute forward and backward scattered field
Ep = Fp*I;
Em = Fm*I;

%% V: modal analysis

% tbd

%% VI: post processing and plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% construct extended mesh (solely for plotting)
pmesh = onedim.generatePlottingMesh(mesh,k0zext/k0,Next);

% compute incident and scattered fields everywhere (*p quantities)
Eip = exp(-1j*k0z*pmesh.z);
EsR = Ep*exp(-1j*k0z*pmesh.z(pmesh.z>mesh.h));
EsL = Em*exp(1j*k0z*pmesh.z(pmesh.z<0));
EsC = -(Z0*I)./mesh.dz;
Esp = [EsL, EsC.', EsR];

% compute total fields everywhere
Etp = Eip+Esp;

% package up currents in extended mesh format
Ip = [EsL*NaN,I.',EsR*NaN];


%Clac Validation case
 results = onedim.slabValidationCase(e0,thetai,k0,k0h,pmesh.z);

%% Error
 results.J(isnan(results.J))=0;
 Ip(isnan(Ip))=0;
 %results.J = [results.J, zeros(1, length(Ip) - length(results.J))];
 Error = sum(abs((Ip - results.J)).^2)/sum(abs(results.J).^2);
 Error_(index) = Error;

end


loglog(N_, Error_, 'o--')
hold on
xlabel('Number of unknowns within slab')
ylabel('Error')
leg = legend ('3-1j*0.2','6-1j*0.2', '5', '1j*0.20', '-1j*100')
[leg,att] = legend('show');
title(leg,'e0 values')
leg.Title.Visible = 'on';


end

