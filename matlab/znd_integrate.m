function [dx, dt] = znd_integrate(lambda,x,P_0, rho_0, T_0, Q, gamma,Rsp,c_0,D,Ea,k)
  [P, U, c, rho,V,x,M,T] = overdrivenD(P_0, rho_0, T_0, Q, gamma,Rsp,c_0,D,lambda);
  dx = (U.-D)/((1-lambda).*k.*e.^(-Ea*T_0./(T)));
  dt = 1.0/((1-lambda).*k.*e.^(-Ea*T_0./(T)));
endfunction

