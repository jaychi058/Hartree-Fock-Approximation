% Eq. (A.11)

function T_pq = T_kinetic(alpha, beta, Rab2)

    ab = alpha * beta;  % alhpa * beta
    a_b = alpha + beta; % alpha + beta
    
    T_pq = ab/(a_b) * (3.0-2.0*ab/a_b*Rab2) * (pi/a_b)^1.5 * exp(-ab/a_b*Rab2);
end