function [D_cj,P_cj,rho_cj, V_cj, C_cj, M_cj, T_cj, u_cj] = znd_polyfn(P_0, rho_0, T_0, Q, gamma,Rsp,c0)
%function evaluates CJ properties
  #M_cj = (2 * q/R*T_0)
  phi = (2 * Q*(((gamma^2)-1)/ gamma))^0.5;
  #p(x) = m^2 - 1 - phi * m.
  p = [1, -phi, -1];
  a = roots (p);
  M_cj  = max(a(1),a(2));
  P_cj = (1 + (gamma * (M_cj^2)))/(1 + gamma);
  P_cj = P_cj * P_0;
  V_cj = (1 + (gamma * (M_cj ^2)))/ (M_cj^2 *( 1 + gamma));
  V_cj = V_cj/rho_0;
  T_cj = (P_cj * V_cj)/Rsp;
  rho_cj = 1/V_cj;
  C_cj = sqrt(gamma * Rsp*1000 * T_cj);
  D_cj = M_cj * c0;
  u_cj = D_cj*(1-rho_0/rho_cj);

endfunction
