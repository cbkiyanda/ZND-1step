function [P, U, c, rho,V,x,M,T] = overdrivenD(P_0, rho_0, T_0, Q, gamma,Rsp,c0,D,lambda)
  %function evalutes parameters of an overdriven function.
  M = D/c0;
  x = (((M - M.^-1)^2 - 2.*(gamma.^2 - 1) .* (lambda .* Q / gamma))./((gamma + 1).^2 .* M.^2)).^0.5;
  V = ((gamma + (M^-2))/(gamma + 1)) - x;
  V = V/rho_0;
  P = (((gamma*M^2) + 1)/(gamma + 1)) + gamma * M^2 * x;
  P = P * P_0;
  U = ((M - (M^-1))/(gamma + 1)) + M*x;
  U = U* c0;
  rho = 1./V;
  T = P .* V./Rsp;
  c = sqrt(gamma.*Rsp.*1090.*T);
endfunction
