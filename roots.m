function F = roots(x, REF)

%A, alpha, delta, rho, theta

F(1) = ((REF(1) * REF(2) *  x(2)^(REF(2) - 1) - REF(4)) / ...
    REF(5)) * x(1); % c dot

F(2) = REF(1) * x(2)^REF(2) - REF(3) * x(2) - x(1); % k dot
end

