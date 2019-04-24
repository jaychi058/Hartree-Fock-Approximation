% Comparsion of the least-square fit of a 1s Slater function obtain at
% STO-1G, STO-2G, and STO-3G level

function [alhpa_new_1, alhpa_new_2, d1, d2] = basisFun(L, zeta1, zeta2, plotFig)

    n_ptr = 101;        % No. of grid points
    Ra = linspace(0, 4, n_ptr);
    r0 = 0.0;
    abs_r_Ra = abs(r0-Ra);
    r_Ra2 = abs_r_Ra .* abs_r_Ra;

    % Coefficients for basis function
    % Eq. (3.219 to 3.221)
    d_coeff = [1.0 0.0 0.0; 0.678914 0.4303129 0.0; 0.444635 0.535328 0.154329];

    % alpha: the Gaussian orbital exponent
    alpha = [0.270950 0.0 0.0; 0.151623 0.851819 0.0; 0.109818 0.405771 2.227660];

    % when zeta = 1.24, the results for alpha_new_2 is as same as Eq. (3.225)
    alhpa_new_1 = (zeta1*zeta1) * alpha;    % Eq. (3.224)
    alhpa_new_2 = (zeta2*zeta2) * alpha;    % Eq. (3.224)

    % Eq. (3.202)
    phi_1s_SF = sqrt(zeta2^3/pi) .* exp(-zeta2*abs_r_Ra);

    % Eq. (3.203)
    d1 = zeros(L, L);
    d2 = zeros(L, L);
    for p = 1:L
        d1(1, p) = d_coeff(1, p) * (2.0*alhpa_new_1(1, p)/pi)^0.75;
        d1(2, p) = d_coeff(2, p) * (2.0*alhpa_new_1(2, p)/pi)^0.75;
        d1(3, p) = d_coeff(3, p) * (2.0*alhpa_new_1(3, p)/pi)^0.75;

        d2(1, p) = d_coeff(1, p) * (2.0*alhpa_new_2(1, p)/pi)^0.75;
        d2(2, p) = d_coeff(2, p) * (2.0*alhpa_new_2(2, p)/pi)^0.75;
        d2(3, p) = d_coeff(3, p) * (2.0*alhpa_new_2(3, p)/pi)^0.75;
    end
    
    % Eq. (3.212, 3.213 to 3.215)
    phi_STO1G_CGF = zeros(1, n_ptr);
    phi_STO2G_CGF = zeros(1, n_ptr);
    phi_STO3G_CGF = zeros(1, n_ptr);
    for p = 1:L
        phi_STO1G_CGF = phi_STO1G_CGF + d2(1, p) .* exp(-alhpa_new_2(1, p)*r_Ra2);
        phi_STO2G_CGF = phi_STO2G_CGF + d2(2, p) .* exp(-alhpa_new_2(2, p)*r_Ra2);
        phi_STO3G_CGF = phi_STO3G_CGF + d2(3, p) .* exp(-alhpa_new_2(3, p)*r_Ra2);
    end

    if(plotFig == true)
        % The result is as same as Fig. 3.3 on page 158 if zeta = 1.0 
        plot(Ra, phi_1s_SF, '-', Ra, phi_STO1G_CGF, ':',...
            Ra, phi_STO2G_CGF, '-.', Ra, phi_STO3G_CGF, '--')
        legend({'SLATER', 'STO-1G', 'STO-2G', 'STO-3G'}, 'FontSize', 16)
        xlabel('Radius (a.u.)', 'FontSize', 16)
        ylabel('\phi_{1s}', 'FontSize', 16)
    end
        
end
