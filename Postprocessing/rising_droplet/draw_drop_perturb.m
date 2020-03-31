%draw droplet and perturbation

close all
clear variables

set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',2,'defaultpatchlinewidth',.7);

N = 1e5;    %numkber of point
theta = 0:pi/(N-1):pi;
modes = 10;        %how many modes
R0 = ones(1,numel(theta));  %unpertirbed interface
%delta = 0.1;  %amplitude of the perturbation

plotINT = 1;

%manyDelta = logspace(1e-8,1,8);
manyDelta = [1e-4 1e-3 1e-2 1e-1 1e0];
errV = zeros(1,numel(manyDelta));
errCons = zeros(1,numel(manyDelta));
Vone0 = zeros(1,numel(manyDelta));
Vone1 = zeros(1,numel(manyDelta));
Vone2 = zeros(1,numel(manyDelta));
Vone3 = zeros(1,numel(manyDelta));

%integration
INT = ([diff(theta),0]+[0,diff(theta)])/2;

%starting volume
V1 = 2/3*pi*INT*(R0.^3.*sin(theta))';
Van = 4/3*pi;
err1 = (V1-Van)/Van;
%display(err1)

for k = 1:numel(manyDelta)
    
    display(k)

    R1 = R0;
    legSUM = 0;
    for i = 2:modes+1
        %legendre function
        PPP = legendre(i,cos(theta));
        %legendre polynomia
        P = PPP(1,:);
        
%         figure(2)
%         hold on
%         plot(P)
%         hold on

        %mode coefficient: I impose it respecting convergence
        f = 1/i^4;

        delta = manyDelta(k);
        R1 = R1 + delta*(2*i+1)*f*P;
        
        legSUM = legSUM + delta*(2*i+1)*f*P;

    end
    
    if plotINT==1
    figure(1)
    subplot(3,2,k)
    plot(R1.*sin(theta),R1.*cos(theta),'r')
    hold on
    plot(-R1.*sin(theta),R1.*cos(theta),'r')
    plot(R0.*sin(theta),R0.*cos(theta),'k--')
    plot(-R0.*sin(theta),R0.*cos(theta),'k--')
    
    axis equal
    axis([-1.5 1.5 -1.5 1.5])
    xlabel('r')
    ylabel('x')
    hold off
    end
    
    %compute volume
    V = 2/3*pi*INT*(R1.^3.*sin(theta))';
    errV(k) = (V-Van)/Van;
    
    %conserved volume
    %Vcons = 2/3*pi*INT*(1 + 3*legSUM + legSUM.^3 + 3*legSUM.^2)';
    Vone0(k) = 2/3*pi*INT*(ones(1,numel(legSUM)).*sin(theta))';
    Vone1(k) = 2/3*pi*INT*(3*legSUM.*sin(theta))';
    Vone2(k) = 2/3*pi*INT*(3*legSUM.^2.*sin(theta))';
    Vone3(k) = 2/3*pi*INT*(legSUM.^3.*sin(theta))';
    %Vcons = 2/3*pi*INT*(1*sin(theta) + 3*legSUM.*sin(theta) +
    %legSUM.^3.*sin(theta))'; DON NOT CONSERVE BECAUSE THERE ARE
    %INTERACTIONS OF QUADRATIC TERMS??!!!
    %errCons(k) = (Vcons-Van)/Van;

end

figure
loglog(manyDelta,abs(errV),'o-')
xlabel('\delta')
ylabel('err_V')
grid on

% figure
% loglog(manyDelta,abs(errCons),'o-')
% xlabel('\delta')
% ylabel('V conserv')
% grid on

% figure
% loglog(manyDelta,Vone0,'o-')
% xlabel('\delta')
% ylabel('V0')
% grid on
% 
% figure
% loglog(manyDelta,Vone1,'o-')
% xlabel('\delta')
% ylabel('V1')
% grid on
% 
% figure
% loglog(manyDelta,Vone2,'o-')
% xlabel('\delta')
% ylabel('V2')
% grid on
% 
% figure
% loglog(manyDelta,Vone3,'o-')
% xlabel('\delta')
% ylabel('V3')
% grid on










