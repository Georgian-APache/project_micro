function h = rician_channel(N, K_dB)
% RICIAN_CHANNEL  Generate N i.i.d. Rician fading coefficients.
%
% Inputs:
%   N    — number of channel coefficients to generate
%   K_dB — Rician K-factor in dB (ratio of LOS power to scattered power)
%
% Output:
%   h — N x 1 complex vector of unit-power Rician samples.
%
% Model: h = sqrt(K/(K+1)) * exp(j*phi_LOS) + sqrt(1/(K+1)) * CN(0,1)
%   where phi_LOS is a uniform random phase per coefficient.
% Total expected power E[|h|^2] = 1.

    K   = 10^(K_dB/10);                                  % linear K-factor
    los = sqrt(K/(K+1)) * exp(1j*2*pi*rand(N,1));        % LOS component
    nlos = sqrt(1/(2*(K+1))) * (randn(N,1) + 1j*randn(N,1));  % scattered
    h   = los + nlos;
end
