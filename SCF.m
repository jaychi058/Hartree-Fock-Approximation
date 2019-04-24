function [E0, E_tot, eplson, C, P] = SCF(N, R, Za, Zb, P, TE, H, X, print)

    % Convergence critera
    delta = 99999;
    crit = 1e-10;
    
    % Maximum number of iterations
    Maxit = 2000; 
    Iter = 0;
    
    E0 = 0.0;           % total electronic energy
    E_tot = 0.0;        % total energy
    
    % ********** Step 5 ********** %
    if (print==true)
        fprintf('Iter.\tP11\tP12\tP22\t   E0(a.u.)\tDelta(Norm2(P-P_old))\n');
    end
    while(crit < delta)
        
        G = zeros(N);
        
        for mu = 1:N
            for nu = 1:N
                for lambda = 1:N
                   for sigma = 1:N
                       G(mu, nu) = G(mu, nu) + P(lambda, sigma)*(TE(mu, nu, sigma, lambda)...
                           -0.5*TE(mu, lambda, sigma, nu));
                   end
                end
            end
        end
        
        % ********** Step 6 ********** %
        F = H + G;      % Eq. (3.154)
        
        E0 = 0.5 * sum(sum(P' .* (H+F)));   % Eq. (3.184)
        
        % ********** Step 7 ********** %
        F_new = X' * F * X;                 % Eq. (3.177)
        
        % ********** Step 8 ********** %
        [C_new, eplson] = eig(F_new);       % Eq. (3.178)
        
        % ********** Step 9 ********** %
        C = X * C_new;                      % Eq. (3.174)
        
        % ********** Step 10 ********** %
        % Eq. (3.145)
        P_old = P;
        P = zeros(N);
        for mu = 1:N
            for nu = 1:N
               for a = 1:(N/2)
                  P(mu, nu) = P(mu, nu) + 2.0*C(mu, a)*C(nu, a);
               end
            end
        end
        
        % ********** Step 11 ********** %
        delta = norm(P - P_old);
        
        Iter = Iter + 1;
        
        if (print == true)
            fprintf('%d\t%6.4f\t%6.4f\t%6.4f\t%12.6f\t%13.11f\n', Iter, ...
                P(1, 1), P(1, 2), P(2, 2), E0, delta);
        end
        % If Iter > Maxit; stop SCF
        if (Iter > Maxit)
            fprintf('\n!!!!!!!!!! The result might not be correct!!!!!!!!!!\n');
            break;
        end
        
    end
    
    % total Energy
    E_tot = E0 + (Za*Zb)/R;             % Eq. (3.185)
    
end