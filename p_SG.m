function p = p_SG(gamma,rho, e, p0)
% % % % EOS SPECIFIC FUNCTION: STIFFENED GAS % % % % %
% Output: pressure (p)
% Inputs: specific heat ratio (gamma), denisity (rho), internal energy per unit mass (e),
% pressure constant (p0->p0 = 0 for ideal gas EOS)
p = e .* (rho .* (gamma-1)) - gamma .* p0;
end