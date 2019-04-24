% Eq. (A.33)

function V_pq = V_nuclear(alpha, beta, Zc, Rab2, Rpc2)

    ab = alpha * beta;  % alhpa * beta
    a_b = alpha + beta; % alpha + beta
    
    F0 = errFunction(a_b*Rpc2);
    
    V_pq = -2.0*pi/a_b * Zc * exp(-ab/a_b*Rab2) * F0;
end