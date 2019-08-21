function [dL, dt] = znd_integrate_x(x,y,P_0, rho_0, T_0, Q, gamma,Rsp,c_0,D,Ea,k)
%function evaluates reaction progress as a function of space.
  [P, U, c, rho,V,x,M,T] = overdrivenD(P_0, rho_0, T_0, Q, gamma,Rsp,c_0,D,y(1));
  dL = ((1-y(1)).*k.*e.^(-Ea*T_0./(T)))/(U.-D);
  dt = 1.0/(U.-D);
  %fprintf('%g %g %g %g %g %g %g %g %g %g \n',P,U,c,rho,V,x,M,T,dL,dt)
endfunction
