
function [S_overlap, H_core, TE, E_H, E_He_p, E_He] = molInteg(N, L, alhpa_new_1, alhpa_new_2, d1, d2, R, Za, Zb)

    S_overlap = eye(N);
    T_kinetic_Mat = zeros(N);
    V1_coulomb_Mat = zeros(N);
    V2_coulomb_Mat = zeros(N);
    H_core = zeros(N);

    TE = zeros(N, N, N, N);

    r0 = 0.0;
    R2 = (r0-R)*(r0-R);
    
    % Solving for the overlap matrix S
    % Eq. (3.228) for S_uv
    % Solving for the kinetic energy matrix T
    % Eq. (3.151) for T_pq
    % Solving for the nuclear attraction matrix V
    % Eq. (3.152) for V_pq
    for p = 1: L
        for q = 1: L
            % the overlap matrix S 
            S_overlap(1, 2) = S_overlap(1, 2) + S_overlap_Mat(alhpa_new_1(L, p), alhpa_new_2(L, q), R2) * d1(L, p) * d2(L, q);

            % the kinetic energy matrix T
            T_kinetic_Mat(1, 1) = T_kinetic_Mat(1, 1) + T_kinetic(alhpa_new_1(L, p), alhpa_new_1(L, q), 0.0) * d1(L, p) * d1(L, q);
            T_kinetic_Mat(1, 2) = T_kinetic_Mat(1, 2) + T_kinetic(alhpa_new_1(L, p), alhpa_new_2(L, q), R2) * d1(L, p) * d2(L, q);
            T_kinetic_Mat(2, 2) = T_kinetic_Mat(2, 2) + T_kinetic(alhpa_new_2(L, p), alhpa_new_2(L, q), 0.0) * d2(L, p) * d2(L, q);

            % the nuclear attraction matrix V
            % Eq. (3.210)
            Rp = (alhpa_new_2(L, p)*r0 + alhpa_new_2(L, q)*R) / (alhpa_new_1(L, p)+alhpa_new_2(L, q));
            Rp2 = Rp * Rp;
            Rc = R - Rp;
            Rc2 = Rc * Rc;

            V1_coulomb_Mat(1, 1) = V1_coulomb_Mat(1, 1) + V_nuclear(alhpa_new_1(L, p), alhpa_new_1(L, q), Za, 0.0, 0.0) * d1(L, p) * d1(L, q);
            V1_coulomb_Mat(1, 2) = V1_coulomb_Mat(1, 2) + V_nuclear(alhpa_new_1(L, p), alhpa_new_2(L, q), Za, R2, Rp2) * d1(L, p) * d2(L, q);
            V1_coulomb_Mat(2, 2) = V1_coulomb_Mat(2, 2) + V_nuclear(alhpa_new_2(L, p), alhpa_new_2(L, q), Za, 0.0, R2) * d2(L, p) * d2(L, q);

            V2_coulomb_Mat(1, 1) = V2_coulomb_Mat(1, 1) + V_nuclear(alhpa_new_1(L, p), alhpa_new_1(L, q), Zb, 0.0, R2) * d1(L, p) * d1(L, q);
            V2_coulomb_Mat(1, 2) = V2_coulomb_Mat(1, 2) + V_nuclear(alhpa_new_1(L, p), alhpa_new_2(L, q), Zb, R2, Rc2) * d1(L, p) * d2(L, q);
            V2_coulomb_Mat(2, 2) = V2_coulomb_Mat(2, 2) + V_nuclear(alhpa_new_2(L, p), alhpa_new_2(L, q), Zb, 0.0, 0.0) * d2(L, p) * d2(L, q);

        end 
    end
    
    S_overlap(2, 1) = S_overlap(1, 2);
    %S_overlap;

    T_kinetic_Mat(2, 1) = T_kinetic_Mat(1, 2);
    T_kinetic_Mat;

    V1_coulomb_Mat(2, 1) = V1_coulomb_Mat(1, 2);
    V1_coulomb_Mat;

    V2_coulomb_Mat(2, 1) = V2_coulomb_Mat(1, 2);
    V2_coulomb_Mat;
    
    H_core = T_kinetic_Mat + V1_coulomb_Mat + V2_coulomb_Mat;

    % Solving for two-electron integrals
    % Eq. (3.155) for two-electron integrals
    for p = 1 : L
        for q = 1 : L
            for r = 1 : L
                for s = 1 : L
                    Rp = alhpa_new_2(L, r)*R / (alhpa_new_2(L, r)+alhpa_new_1(L, s));
                    Rp2 = Rp * Rp;
                    Rt = alhpa_new_2(L, p)*R / (alhpa_new_2(L, p)+alhpa_new_1(L, q));
                    Rt2 = Rt * Rt;

                    Rpq = R - Rp;
                    Rpq2 = Rpq * Rpq;

                    Rpt = Rp - Rt;
                    Rpt2 = Rpt * Rpt;

                    TE(1, 1, 1, 1) = TE(1, 1, 1, 1) + Two_E_Integral(alhpa_new_1(L, p), alhpa_new_1(L, q), alhpa_new_1(L, r), alhpa_new_1(L, s), 0.0, 0.0, 0.0) ...
                        * d1(L, p) * d1(L, q) * d1(L, r) * d1(L, s);
                    TE(2, 2, 1, 1) = TE(2, 2, 1, 1) + Two_E_Integral(alhpa_new_2(L, p), alhpa_new_2(L, q), alhpa_new_1(L, r), alhpa_new_1(L, s), 0.0, 0.0, R2) ...
                        * d2(L, p) * d2(L, q) * d1(L, r) * d1(L, s);
                    TE(2, 1, 1, 1) = TE(2, 1, 1, 1) + Two_E_Integral(alhpa_new_2(L, p), alhpa_new_1(L, q), alhpa_new_1(L, r), alhpa_new_1(L, s), R2, 0.0, Rt*Rt) ...
                        * d2(L, p) * d1(L, q) * d1(L, r) * d1(L, s);
                    TE(2, 2, 2, 1) = TE(2, 2, 2, 1) + Two_E_Integral(alhpa_new_2(L, p), alhpa_new_2(L, q), alhpa_new_2(L, r), alhpa_new_1(L, s), 0.0, R2, Rpq2) ...
                        * d2(L, p) * d2(L, q) * d2(L, r) * d1(L, s);
                    TE(2, 1, 2, 1) = TE(2, 1, 2, 1) + Two_E_Integral(alhpa_new_2(L, p), alhpa_new_1(L, q), alhpa_new_2(L, r), alhpa_new_1(L, s), R2, R2, Rpt2) ...
                        * d2(L, p) * d1(L, q) * d2(L, r) * d1(L, s);
                    TE(2, 2, 2, 2) = TE(2, 2, 2, 2) + Two_E_Integral(alhpa_new_2(L, p), alhpa_new_2(L, q), alhpa_new_2(L, r), alhpa_new_2(L, s), 0.0, 0.0, 0.0) ...
                        * d2(L, p) * d2(L, q) * d2(L, r) * d2(L, s);
                end
            end
        end
    end
    
    TE(1, 1, 1, 1);

    TE(2, 2, 1, 1);
    TE(1, 1, 2 ,2) = TE(2, 2, 1, 1);

    TE(2, 1, 1, 1);
    TE(1, 2, 1, 1) = TE(2, 1, 1, 1);
    TE(1, 1, 2, 1) = TE(2, 1, 1, 1);
    TE(1, 1, 1, 2) = TE(2, 1, 1, 1);

    TE(2, 2, 2, 1);
    TE(2, 2, 1, 2) = TE(2, 2, 2, 1);
    TE(2, 1, 2, 2) = TE(2, 2, 2, 1);
    TE(1, 2, 2, 2) = TE(2, 2, 2, 1);

    TE(2, 1, 2, 1);
    TE(1, 2, 2, 1) = TE(2, 1, 2, 1);
    TE(2, 1, 1, 2) = TE(2, 1, 2, 1);
    TE(1, 2, 1 ,2) = TE(2, 1, 2, 1);

    TE(2, 2, 2, 2);
    E_H = T_kinetic_Mat(2,2)+V2_coulomb_Mat(2,2);
    E_He_p = T_kinetic_Mat(1,1)+V1_coulomb_Mat(1,1);
    E_He = 2.0*(T_kinetic_Mat(1,1)+V1_coulomb_Mat(1,1))+TE(1, 1, 1, 1);

end