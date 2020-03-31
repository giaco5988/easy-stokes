%compute non dimensional parameters for experiments

%fluid 1 is th einner fluid, fluid 2 is the outer fluid
%give the physical parameters
mu1 = 1e-3;                 % inner viscosity
mu2 = 39.1;                 % outer viscosity
rho1 = 800;                   % inner density
rho2 = 1021;                % outer density
g = 9.8;                    % gravity
gamma = 30*1e-3;            % surface tension
Rin = 5e-3;                 % radius as an input
%CaIn = 1;                  % capillary number as an input

%settling velocirties (Batchelor pag. 236)
V = -1/3*(Rin^2*g/mu2)*(rho1-rho2)*(mu2+mu1)/(mu2+1.5*mu1);
display(V)

%compute capillary number as Leal
Ca = mu2*V/gamma;
display(Ca)

%compute radius form fiven capillary number
