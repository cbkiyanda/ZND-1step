function [dx, dt] = znd_integrate(lambda,xt,P_0, rho_0, T_0, Q, gamma,Rsp,c_0,D,Ea,k)
  [P, U, c, rho,V,x,M,T] = overdrivenD(P_0, rho_0, T_0, Q, gamma,Rsp,c_0,D,lambda);
  dx = (U.-D)/((1.0-lambda).*k.*e.^(-Ea*T_0./(T)));
  dt = dx/(U.-D);
  %fprintf('%g %g %g %g %g %g %g %g %g %g %g \n',P,U,c,rho,V,x,M,T,U-D,dx,dt)
endfunction

