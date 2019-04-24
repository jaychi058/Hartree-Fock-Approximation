function F0 = errFunction(t)
        if (t < 1e-6)
            F0 = 1.0 - t/3.0;
        else
            F0 = 0.5 * sqrt(pi/t) * erf(sqrt(t));
        end
end