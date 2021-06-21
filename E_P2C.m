function E = E_P2C(rho, u, e)
% Performing a primitive to conservative operation
% Evaluate the total energy per unit volume (E)
% Inputs: internal energy per unit mass (e), density (rho), velocity (u)
E = (e + (u.^2)/2).*rho;
end