% Eq. (A.41)

function ABCD = Two_E_Integral(alpha, beta, gamma, delta, Rab2, Rcd2, Rpq2)
    ab = alpha * beta;  % alhpa * beta
    a_b = alpha + beta; % alpha + beta
    rd = gamma * delta;     % gamma * delta
    r_d = gamma + delta;    % gamma + delta
    
    F0 = errFunction(a_b*r_d/(a_b+r_d)*Rpq2);
    
    ABCD = 2.0*pi^2.5 / ((a_b)*(r_d)*sqrt(a_b+r_d)) ...
        * exp(-ab/a_b*Rab2 - rd/r_d*Rcd2) * F0;

end