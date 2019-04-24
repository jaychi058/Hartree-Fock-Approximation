clc; clear; close all;

% ********** The Hartree-Fock SCF procedure is on the page 146 ********** %
% ********** 3.5.2 STO-3G H2 **********

% ********** Step 1 ********** %
% zeta2 = 1.0, the plot result of phi is as same as Fig. (3.3) on page 158
R = 1.40;               % <<<<<========== change R here ==========
Za = 1.0;               % Za = 1: H
Zb = Za;                % Zb = 1: H
zeta2 = 1.24;           % <<<<<========== change zeta here ==========
zeta1 = zeta2;

N = 2;                  % Number of electorns
L = 3;                  % the Length of linear combination for STO
                        % 1: STO-1G
                        % 2: STO-2G
                        % 3: STO-3G

STO = sprintf('STO-%dG', L);

[alhpa_new_1, alhpa_new_2, d1, d2] = basisFun(L, zeta1, zeta2, true);

% ********** Step 2 ********** %
[S_uv, H_core_uv, TE, E_H, E_He_p, E_He] = molInteg(N, L, alhpa_new_1, alhpa_new_2, d1, d2, R, Za, Zb);

% ********** Step 3 ********** %
X_mat = canonicalOrth(N, S_uv);

% ********** Step 4 ********** %
P_guess = rand(N);      % Initial Guess density matrix

% ********** After Step 5 (SCF Loop) ********** %
[E0, E_tot, eplson, C, P] = SCF(N, R, Za, Zb, P_guess, TE, H_core_uv, X_mat, true);

fprintf('\nThe final converged total Electronic Energy(a.u.): %13.10f\n', E0);
fprintf('\nThe H atom Energy(a.u.): %13.10f\n', E_H);
fprintf('\nThe final converged total Energy(a.u.): %13.10f\n', E_tot);

% ********** double check N by using Mulliken population ********** %
fprintf('\nUsing N =Tr(PS) to double check N:')
N_check = trace(P * S_uv)% Eq. (3.195)


% ********************************************************************* %
% **********     This part is for the RHF potential curve    ********** %
% ********** The rsult is as smae as the fig. 3.5 (Page 166) ********** %
% ********************************************************************* %
R = 0.1:0.1:4.2;
lenR = numel(R);
E_H = 0.0;
E_tot = zeros(1, lenR);

[alhpa_new_1, alhpa_new_2, d1, d2] = basisFun(L, zeta1, zeta2, false);

% ********** Step 2 ********** %
for i = 1:lenR
    [S_uv, H_core_uv, TE, E_H, E_He_p, E_He] = molInteg(N, L, alhpa_new_1, alhpa_new_2, d1, d2, R(i), Za, Zb);

    % ********** Step 3 ********** %
    X_mat = canonicalOrth(N, S_uv);

    % ********** Step 4 ********** %
    P_guess = rand(N);      % Initial Guess density matrix

    % ********** After Step 5 (SCF Loop) ********** %
    [E0, E_tot(i), eplson, C, P] = SCF(N, R(i), Za, Zb, P_guess, TE, H_core_uv, X_mat, false);
end

figure(2)
plot(R, E_tot-2.0*E_H)
xlim([0.0, 4.2])
ylim([-0.3, 0.6])
xlabel('Radius (a.u.)', 'FontSize', 16)
ylabel('E(H_{2}) - 2E(H) (a.u.)', 'FontSize', 16)
legend({STO}, 'FontSize', 16)
grid on