function dY_dt = ode(t, ode_initial, REF)

% A, alpha, delta, rho, theta

c = ode_initial(1, 1);
k = ode_initial(2, 1);

dY_dt(1, 1) = (REF(1) * REF(2) *  k^(REF(2) - 1) - ...
    REF(3) - REF(4)) / REF(5) * c; % dc/dt

dY_dt(2, 1) = REF(1) * k^REF(2) - REF(3) * k - c; % dk/dt
end

% t will be used in ode45