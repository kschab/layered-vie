close all
clear
clc

% example: homogeneous slab layer validation
% 
% this script simulates a homogeneous dielectric layer and compares results
% to analytic models
%
% scu em lab
% 2021

%% I: user defined parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% excitation parameters
thetai = pi/3;
k0 =  1;
e0 = 3-1j*0.2;     % just a parameter defining arbitrary permittivity profiles

% geometry parameters
k0h = 6.4;  % free space electrical thickness of slab
N = 100;    % number of unknowns within slab
k0zext = [-k0h, 2*k0h]; % free space electrical size of extended mesh for plotting 
Next = 75; % number of points in each exterior plot region

%% II: setup 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% permittivity profile -- demo of two ways of defining the layered media
epr = ones(N,1)*e0;

% create mesh
h = k0h/k0;
mesh = onedim.generateMesh(N,h);

% assign permittivity
mesh.epr = epr;

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

% plotting fields and currents on extended mesh
figure()
onedim.plot1D(pmesh,real(Eip));
hold on
onedim.plot1D(pmesh,real(Esp));
onedim.plot1D(pmesh,real(Etp));
onedim.plot1D(pmesh,real(Ip));
onedim.geometryOverlay1D(mesh)
xlabel('z, m')
ylabel('a.u.')
title('time domain snapshot')
legend('E_i -- incident field','E_s -- scattered field','E_t -- total field','I -- polarization current','location','best')

% plotting fields and currents on extended mesh
figure()
subplot(2,1,1)
a = onedim.plot1D(pmesh,abs(Eip));
hold on
b = onedim.plot1D(pmesh,abs(Esp));
c = onedim.plot1D(pmesh,abs(Etp));
onedim.geometryOverlay1D(mesh);
xlabel('z, m')
ylabel('|E|, a.u.')

%% VII: analytic validation case
results = onedim.slabValidationCase(e0,thetai,k0,k0h,pmesh.z);

subplot(2,1,1)
d = plot(results.kz/k0,abs(results.Ei),'o');
f = plot(results.kz/k0,abs(results.Es),'o');
g = plot(results.kz/k0,abs(results.Et),'o');
legend([a,b,c,d,f,g],{'E_i -- incident field (MoM)','E_s -- scattered field (MoM)','E_t -- total field (MoM)','E_i -- incident field (Analytic)','E_s -- scattered field (Analytic)','E_t -- total field (Analytic)'},'location','best');


subplot(2,1,2)
onedim.plot1D(pmesh,real(Ip));
hold on
plot(results.kz/k0,real(results.J),'o');
onedim.plot1D(pmesh,imag(Ip));
hold on
plot(results.kz/k0,imag(results.J),'o');
ylim([-1.2,1.2]*max(abs(results.J)));
xlabel('z, m')
ylabel('J, a.u.')
legend('Re J (MoM)','Re J (Analytic)','Im J (MoM)','Im J (Analytic)','location','best')