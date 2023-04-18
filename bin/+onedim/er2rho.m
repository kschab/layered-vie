function rho = er2rho(k,er)
%{
er2rho
--------------------
converts complex relative permittivity to resistivity

inputs
    k           free space wavenumber           [scalar, real]
    er          complex relative permittivity   [scalar, complex]

outputs
    rho         complex resistivity             [scalar, complex]

uses the e(jwt) time convention, lossy materials should be represented by
er = e' - je'', with e'' >= 0.  This leads to Re{rho}>=0

%}

constants = onedim.constants;
eta = constants.eta;
rho = 1./(1j*k./eta*(er-1));
end
