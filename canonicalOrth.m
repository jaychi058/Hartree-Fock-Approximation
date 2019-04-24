% I got the Matrix U just from the Eq. (3.259)
% I got s_1_2 just from the Eq. (3.260)
function X_mat = canonicalOrth(N, S_uv)

    X_mat = zeros(N);
    U_mat = zeros(N);
    s_1_2 = zeros(N);
    
    % Eq. (3.259)
    U_mat = [1.0/sqrt(2.0) 1.0/sqrt(2.0); 1.0/sqrt(2.0) -1.0/sqrt(2.0)];
    % Eq. (3.260)
    s_1_2 = [1.0/sqrt(1.0+S_uv(1, 2)) 0.0; 0.0 1.0/sqrt(1.0-S_uv(1, 2))];
    % Eq. (3.262)
    X_mat = U_mat * s_1_2;

end