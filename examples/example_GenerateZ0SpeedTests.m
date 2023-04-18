close all
clear
clc

% example: checks the accuracy and speed of generateZ0Fast
%
% this script examines the relative speed of generateZ0 and an accelerated
% version generateZ0Fast.  the latter is approx. 30x faster for large
% systems.
%
% scu em lab
% 2021


ndex = 0;
N_ = floor(logspace(1,4,11));
t_GenerateZ0 = 0*N_;
t_GenerateZ0Fast = 0*N_;
matrixErrorMax = 0*N_;
matrixErrorAvg = 0*N_;
for N = N_
    ndex = ndex+1;
    %% I: user defined parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % excitation parameters
    thetai = pi/3;
    k0 =  1;
    e0 = 3-1j*0.2;     % just a parameter defining arbitrary permittivity profiles
    
    % geometry parameters
    k0h = 6.4;  % free space electrical thickness of slab
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
    tic()
    Z0 = onedim.GenerateZ0(mesh, k0z, k0);
    t_GenerateZ0(ndex) = toc();
    
    tic()
    Z0F = onedim.GenerateZ0Fast(mesh, k0z, k0);
    t_GenerateZ0Fast(ndex) = toc();
    
    matrixErrorAvg(ndex) = sum(sum(abs(Z0-Z0F)))/N^2;
    matrixErrorMax(ndex) = max(max(abs(Z0-Z0F)));
    
    %% plotting
    figure(100)
    subplot(2,1,1)
    cla()
    loglog(N_,t_GenerateZ0,'o-')
    hold on
    loglog(N_,t_GenerateZ0Fast,'x-')
    loglog(N_,t_GenerateZ0./t_GenerateZ0Fast,'^-')
    grid on
    xlabel('number unknowns')
    ylabel('time, s')
    legend('GenerateZ0','GenerateZ0Fast','speedup','location','best')
    subplot(2,1,2)
    cla()
    loglog(N_,matrixErrorAvg,'o-')
    hold on
    loglog(N_,matrixErrorMax,'x-')
    xlabel('number unknowns')
    ylabel('matrix error')
    legend('mean','max','location','best')
    grid on
    drawnow() 
end

