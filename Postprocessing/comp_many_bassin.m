%plot comparison between linear and non linear simulations at different
%time step

%clear all

set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',1,'defaultpatchlinewidth',.7);

R = 0.06:0.004:0.14;
theta = 0:0.3:2*pi;

a1 = R.*cos(theta)-0.06;
a2 = R.*sin(theta);

stab_unstab = zeros(numel(a1),numel(a2));

parameters = zeros(441,2);
parameters(:,2) = repmat(a2',21,1);
aaa = repmat(a1,21,1);
aaa = aaa(:)';
parameters(:,1) = aaa;

for k = 1:numel(parameters)/2
    
    disp(k)
    
    filename = strcat('q=200_visc=0.5_Dt=0.01_loop=',num2str(end_loop),'_epsilon1=',num2str(parameters(k,1)),'_epsilon2=',num2str(parameters(k,2)),'_Ca=6_RK2.mat');
    load(filename)
    
    time = 1:10:find(risb(2,:)==0,1,'first')-1;
    
    if numel(time)<1
        time = 10:10:end_loop;
        %ccc = num2str(k);
%         if str2double(ccc(end))==0
%             stab_unstab(ceil(k/21),10) = 0;
%         else
            %stab_unstab(ceil(k/21),str2double(ccc(end-1:end))) = 0;
            
            %stab_unstab(ceil(k/21),k-21*(floor((k-0.1)/21))) = 0;
        %end
    else
        stab_unstab(ceil(k/21),k-21*(floor((k-0.1)/21))) = 1;
    end

    elon_mine = zeros(numel(time),1);
    myarea = zeros(numel(time),1);
    V_in = axis_int(risa(:,1)',risb(:,1)');

end

figure
[X,Y] = meshgrid(a1,a2);
surf(X,Y,stab_unstab')
xlabel('a1')
ylabel('a2')
title('bassin of attraction')

figure
[X,Y] = meshgrid(a1,a2);
contour(X,Y,stab_unstab')
xlabel('a1')
ylabel('a2')
title('bassin of attraction')

figure
[X,Y] = meshgrid(a1,a2);
contourf(X,Y,stab_unstab')
xlabel('a1')
ylabel('a2')
title('bassin of attraction')