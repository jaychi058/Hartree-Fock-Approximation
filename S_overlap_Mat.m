% Eq. (A.9)

function S_pq = S_overlap_Mat(alpha, beta, Rab2)

    S_pq = (pi./(alpha+beta)).^1.5 .* exp(-alpha.*beta./(alpha+beta)*Rab2);
    
end