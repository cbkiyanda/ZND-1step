%Chemical parameters
gamma = 1.2;
Q   = 50; %Non-dimensionalized
Ea  = 50; %Non-dimensionalized
M = 28.96; #Kg/KMol

%Overdrive = (D/D_cj)^2
f = 1.6; 

%Quiescent State
P_0 = 100;  # kPa
T_0 = 298;  # K

%Constants
k = 1;     % 1/s, Set k to unity for initial run
R = 8.314; %kJ/kMol.K
Rsp   = R / M; # kJ/kg.K
rho_0 = P_0 / (Rsp *T_0) #kg/m^3;
c_0 = sqrt(gamma * Rsp*1000 * T_0);

%
lambda = 1; % End of reaction zone
N12 = 100;  % Number of points per half-reaction zone length

%Dimensionalize heat release and activation energy
q = Q* Rsp *1000*T_0; # KJ/Kg
ea= Ea* Rsp *1000*T_0; # KJ/Kg


[D_cj,P_cj,rho_cj, V_cj, C_cj, M_cj, T_cj, u_cj] = znd_polyfn(P_0, rho_0, T_0, Q, gamma,Rsp,c_0);
D = sqrt(f)*D_cj;

  tspan = [0 0.5 0.9999]
  opt = odeset ("RelTol", 1e-12, "InitialStep", 1e-1, "MaxStep", 1e-1);#, "AbsTol", 1e-8);
  [lambda, xt] = ode45 (@(lambda,xt) znd_integrate(lambda,xt,P_0, rho_0, T_0, Q, gamma,Rsp,c_0,D,Ea,k), tspan, [0,0]);
  t12 = xt(2,2)
  x12 = xt(2,1)
  xr  = xt(3,1)
  
  dx  = x12/N12;
  
  xspan = 0:dx:xr;
  #N   = floor(xr/dx)+1; 
  #xspan = linspace(0,xr,N);
  opt = odeset ("RelTol", 1e-12, "InitialStep", 1e-1, "MaxStep", 1e-1);#, "AbsTol", 1e-8);
  [xsol, lt] = ode45 (@(x,lt) znd_integrate_x(x,lt,P_0, rho_0, T_0, Q, gamma,Rsp,c_0,D,Ea,k), xspan, [0,0]);
  x = xsol;
  lambda = lt(:,1);
  t      = lt(:,2);
  
  P   = 0*lambda;
  U   = 0*lambda;
  c   = 0*lambda;
  rho = 0*lambda;
  V   = 0*lambda;
  xFD = 0*lambda;
  M   = 0*lambda;
  T   = 0*lambda;
  j   = length(lambda);

for i = 1:j; 
  [P(i), U(i), c(i), rho(i),V(i),xFD(i),M(i),T(i)] = overdrivenD(P_0, rho_0, T_0, Q, gamma,Rsp,c_0,D,lambda(i));
end
  rate = k*(1-lambda).*e.^(-Ea*T_0./T);
 figure(1)
 plot (x,P/P_0, "or" )
 grid
 title('Pressure vs x')
 figure(2)
 plot (x, U/c_0)
 grid
 title('Particle Velocity vs x')
 figure(3)
 plot (t/t12,T/T_0, "r" )
 grid
 title('Temperature vs t')
 figure(4)
 plot (x,rho/rho_0, "+bk" )
 grid
 title('Density vs x')
 figure(5)
 plot (x,lambda, "+bk" )
 grid
 title('lambda vs x')
 figure(6)
 plot (t/t12,lambda, "+bk" )
 grid
 title('lambda vs t')
 figure(7)
 plot (x/x12,T/T_0, "r" )
 grid
 title('Temperature vs x')
 figure(8)
 plot(t/t12, rate)
 title('Rate vs time')
  
 # fid_p = fopen("p", "r+t")
 # fid_U = fopen("U", "r+t")
 # fid_Y = fopen("Y", "r+t")
 # fid_T = fopen("T", "r+t")
  
 # fskipl (fid_p,22)
 # fskipl (fid_T,22)
 # fskipl (fid_Y,22)
 # fskipl (fid_U,22)
  
 #  m = length(P);
 #for i = m:-1:1
    %fputs(fid, '%f \n' , P(1:4))
  #  fprintf(fid_p, '%30.15f', P(i)*1000)
    %fputs(fid, '%f \n' , P(1:5))
  #  fskipl (fid_p)
    
  #  fprintf(fid_T, '%30.15f', T(i))
  #  fskipl (fid_T)
    
  #  fprintf(fid_Y, '%30.15f', q*(1-lambda(i)))
  #  fskipl (fid_Y)

   # fprintf(fid_U, '(%30.15f 0 0)', U(i)-U(end))
  #fskipl (fid_U)
    
 #endfor
 # fclose(fid_p)
 # fclose(fid_U)
 # fclose(fid_Y)
 #  fclose(fid_T)  
  

  
%open the files for reading only
%fid_P = fopen('p','rt');
%fid_U = fopen('U','rt');
%fid_T = fopen('T','rt');
%fid_Y = fopen('Y','rt');

%reads the contents of the files into variables. Those variables are
%cell array variables
%P_file_data = textscan(fid_P, '%s', 'Delimiter', '\n');
%U_file_data = textscan(fid_U, '%s', 'Delimiter', '\n');
%T_file_data = textscan(fid_T, '%s', 'Delimiter', '\n');
%Y_file_data = textscan(fid_Y, '%s', 'Delimiter', '\n');

%This is to eliminate one layer of the cell array. I'm not 100% sure why it's
%needed, but it is.
%P_file_data = P_file_data{1};
%U_file_data = U_file_data{1};
%T_file_data = T_file_data{1};
%Y_file_data = Y_file_data{1};


%close the original files
%fclose(fid_P);
%fclose(fid_U);
%fclose(fid_T);
%fclose(fid_Y);

%see how any points are in the znd solution
%m = length(P);

%replace the values with the znd data in the proper lines
%for i = m:-1:1
%  P_file_data{21+m-i} = sprintf('%f', P(i)*1000);
%  U_file_data{21+m-i} = sprintf('(%f 0 0)',U(i)-U(end));
%  T_file_data{21+m-i} = sprintf('%f', T(i));
%  Y_file_data{21+m-i} = sprintf('%f',(((1 - lambda(i))*q)));
%end

%open files for writing. If you use the same filenames, i.e. p, T,Y, U
%the contents of the file will be overwritten
%fid_P = fopen('p', 'wt');
%fid_U = fopen('U', 'wt');
%fid_T = fopen('T', 'wt');
%fid_Y = fopen('Y', 'wt');

%write the file data to the file
%fprintf(fid_P, '%s\n', P_file_data{:});
%fprintf(fid_U, '%s\n', U_file_data{:});
%fprintf(fid_T, '%s\n', T_file_data{:});
%fprintf(fid_Y, '%s\n', Y_file_data{:});

%close the files that you are outputting to
%fclose(fid_P);
%fclose(fid_U);
%fclose(fid_T);
%fclose(fid_Y);

#print -dpng plot1.png
#print -dpng plot1.png


