function [results] = slabValidationCase(epr,thetai,k0,k0h,kz)

    %% set up calculations
    const = onedim.constants();
    chi = epr - 1;  % slab susceptibility
    rhonorm = -imag(chi)/abs(chi)^2 - 1j*real(chi)/abs(chi)^2; % normalized resistivity
    rho = rhonorm/(k0/const.eta);   % full resistivity

    % make sure sign of permittivity square root is chosen correctly
    sepr = sqrt(epr);
    if imag(sepr)>0
        sepr = -sepr;
    end

    % transverse wavenumbers
    kx = k0*sin(thetai);

    % longitudinal wavenumbers
    k0z = sqrt(1-sin(thetai)^2);
    kdz = sqrt(epr-sin(thetai)^2);

    % absolute value of slab thickness
    h = k0h/k0;

    % longitudinal electrical thickness of slab in terms of slab permittivity
    kdzh = kdz*h;
    k0zh = k0z*h;

    % set up z axis
%     z = linspace(min(k0zext),max(k0zext),300)/k0;
%     kz = k0*z;
%kz = k0zext;
z = kz/k0;

    %% system construction and solution
    A = [1, -1, -1, 0;...
        0, exp(-1j*kdzh), exp(1j*kdzh), -exp(-1j*k0zh);...
        k0z, kdz, -kdz, 0;...
        0, kdz*exp(-1j*kdzh), -kdz*exp(1j*kdzh),-k0z*exp(-1j*k0zh)];

    g = [-1; 0 ; k0z; 0];
    f = A\g;
    a = f(1);
    b = f(2);
    c = f(3);
    d = f(4);

    %% reconstruct total fields
    E1 = kz*0;
    E2 = kz*0;
    E3 = kz*0;
    H1 = kz*0;
    H2 = kz*0;
    H3 = kz*0;
    J2 = kz*0;
    for i = 1:length(kz)
        z_ = z(i);
        if z_ <= 0
            E1(i) = exp(-1j*k0z*z_) + a*exp(1j*k0z*z_);
            H1(i) = k0z*(exp(-1j*k0z*z_) - a*exp(1j*k0z*z_));
        elseif z_ <= h
            E2(i) = b*exp(-1j*kdz*z_) + c*exp(1j*kdz*z_);
            H2(i) = kdz*b*exp(-1j*kdz*z_) - kdz*c*exp(1j*kdz*z_);
        else
            E3(i) = d*exp(-1j*k0z*z_);
            H3(i) = k0z*d*exp(-1j*k0z*z_);
        end
    end
    J2 = E2./rho;
    J2(J2==0) = NaN;
    Et = E1+E2+E3;

    %% reconstruct scattered
    Ein = kz*0;
    Esc = kz*0;
    for i = 1:length(kz)
        z_ = z(i);
        Ein(i) = exp(-1j*k0z*z_);
        if z_ <= 0 
            Esc(i) = a*exp(1j*k0z*z_);
        elseif z_ <= h
            Esc(i) = b*exp(-1j*kdz*z_) + c*exp(1j*kdz*z_) - Ein(i);
        else
            Esc(i) = d*exp(-1j*k0z*z_) - Ein(i);
        end
    end

    results.Ei = Ein;
    results.Es = Esc;
    results.Et = Et;
    results.J = J2;
    results.kz = kz;
end