close all
clear
clc

% example: checks the accuracy of the dissimilar media add-on
%
% scu em lab
% 2021



%% I: user defined parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HSER = 30;  % halfspace permittivity

% sweep over incidence angle
tdex = 0;
for thetai = [0, pi/8,pi/4,pi*3/8]
    tdex = tdex+1;
    
    % sweep over frequency, normalized to electrical thickness
    dex = 0;
    K0H = linspace(0.01,1.5,51);
    for k0h = K0H  % free space electrical thickness of slab
        
        dex = dex+1;
        
        e0 = sqrt(HSER);     % slab permittivity according to qwt design rule
        
        k0zext = 2*[-k0h, 2*k0h]; % free space electrical size of extended mesh for plotting
        Next = 75; % number of points in each exterior plot region
        N = 101;
        
        % halfspace parameters
        hs.h = 1;
        hs.er = HSER;
        
        % free space wavenumber
        k0 = k0h/hs.h;
        
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
        
        % halfspace calculations
        eta=376.73;
        eh = hs.er;
        h = hs.h;
        kh = k0*sqrt(eh);
        kzh = sqrt(kh^2-kx^2);
        ct = k0z/k0;
        cth = kzh/(k0*sqrt(eh));
        etah = eta/sqrt(eh);
        
        % halfspace reflection coefficient
        gh = (etah/cth - eta/ct)/(etah/cth + eta/ct);
        
        %% III: calculate impedance matrices, excitation vector, + other operators
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % compute incident field
        V0 = onedim.incidentPlusPlaneWave(mesh,k0z);
        VG = onedim.incidentPlusPlaneWaveGamma(mesh,hs,k0,k0z);
        
        % compute vacuum impedance matrix
        Z0 = onedim.GenerateZ0Fast(mesh, k0z, k0);
        Zp = onedim.GenerateZp(mesh,k0);
        ZG = onedim.GenerateZGamma(mesh, hs, k0z, k0);
        
        % compute forward and backward going wave operators
        [Fm, Fp]= onedim.GenerateFpFmGamma(mesh, hs, k0z, k0);
        
        %% IV: calculate driven solution
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % compute currents
        I = (Z0+Zp+ZG)\(V0+VG);
        
        % compute forward and backward scattered field
        Ep = Fp*I;
        Em = Fm*I;
        
        % reflection and transmission coefficient calculations
        G(dex) = gh*exp(-1j*k0z*2*h)+Em;
        
        %% V: analytic calculation using transmission line approach
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        k1 = k0*sqrt(e0);
        kz1 = sqrt(k1^2 - kx^2);
        
        Zw0 = eta/ct;
        Zw1 = eta/sqrt(e0)*k1/kz1;
        Zw2 = etah/cth;
        
        Zin = Zw1*(Zw2+1j*tan(kz1*h)*Zw1)/(Zw1+1j*tan(kz1*h)*Zw2);
        Ga(dex) = (Zin-Zw0)/(Zin+Zw0);
    end
    
    figure(9)
    subplot(2,1,1)
    cols = get(gca,'colororder');
    plot(K0H,abs(G),'o','color',cols(tdex,:))
    hold on
    plot(K0H,abs(Ga),'color',cols(tdex,:))
    ylim([0,1])
    xlabel('kh')
    ylabel('|\Gamma|')
    legend('MoM','TL model','location','best')    
    subplot(2,1,2)
    cols = get(gca,'colororder');
    plot(K0H,angle(G),'o','color',cols(tdex,:))
    hold on
    plot(K0H,angle(Ga),'color',cols(tdex,:))
    ylim([-pi*1.5,pi*1.5])
    xlabel('kh')
    ylabel('\angle\Gamma')
    legend('MoM','TL model','location','best')
end

%{
%% VI: post processing and plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% construct extended mesh (solely for plotting)
pmesh = onedim.generatePlottingMesh(mesh,k0zext/k0,Next);

% compute incident and scattered fields everywhere (*p quantities)
EiL = exp(-1j*k0z*pmesh.z(pmesh.z<mesh.h))+gh*exp(-1j*2*k0z*h)*exp(1j*k0z*pmesh.z(pmesh.z<mesh.h));
EiR = exp(-1j*k0z*h)*(1+gh)*exp(-1j*kzh*(pmesh.z(pmesh.z>=mesh.h)-h));
EsR = Ep*exp(-1j*k0z*h)*exp(-1j*kzh*(pmesh.z(pmesh.z>mesh.h)-h));
EsL = Em*exp(1j*k0z*pmesh.z(pmesh.z<0));
EsC = -((Z0+ZG)*I)./mesh.dz;
Esp = [EsL, EsC.', EsR];
Eip = [EiL, EiR];

% compute total fields everywhere
Etp = Eip+Esp;

% package up currents in extended mesh format
Ip = [EsL*NaN,I.',EsR*NaN];


% plotting fields and currents on extended mesh
figure()
subplot(2,1,1)
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

subplot(2,1,2)
onedim.plot1D(pmesh,abs(Eip));
hold on
onedim.plot1D(pmesh,abs(Esp));
onedim.plot1D(pmesh,abs(Etp));
xlabel('z, m')
ylabel('a.u.')
title('field magnitudes')
legend('E_i -- incident field','E_s -- scattered field','E_t -- total field','location','best')
%}