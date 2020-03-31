%plot comparison between linear and non linear simulations at different
%time step

%clear all

set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',1,'defaultpatchlinewidth',.7);

visc = 0.5;
Ca = 5;

G = [-0.1 -0.25 -0.5 -1 0.5 0.25 0.1];
%G = -1;
delta = -[0.003:0.001:0.007];
delta = [delta -delta];

stab_unstab = zeros(numel(G),numel(delta));

parameters = zeros(numel(delta)*numel(G),2);
parameters(:,2) = repmat(delta',numel(G),1);
aaa = repmat(G,numel(delta),1);
aaa = aaa(:)';
parameters(:,1) = aaa;

for k = 1:numel(parameters)/2
    
    disp(k)
    
    filename = strcat('G=',num2str(parameters(k,1)),'_q=200_visc=',num2str(visc),'_Dt=0.05_loop=2000_DELTA=',num2str(parameters(k,2)),'_Ca=',num2str(Ca),'_RK2.mat');
    load(filename)
    
    time = 1:10:find(risb(2,:)==0,1,'first')-1;
    
    if numel(time)<1
        time = 10:10:end_loop;
    else
        stab_unstab(ceil(k/numel(delta)),k-numel(G)*(floor((k-0.1)/numel(G)))) = 1;
    end

    elon_mine = zeros(numel(time),1);
    myarea = zeros(numel(time),1);
    V_in = axis_int(risa(:,1)',risb(:,1)');

end

figure
[X,Y] = meshgrid(1:numel(G),delta);
surf(X,Y,stab_unstab')
xlabel('G(t)')
ylabel('\delta^*')

figure
[X,Y] = meshgrid(1:numel(G),delta);
contour(X,Y,stab_unstab')
xlabel('G(t)')
ylabel('\delta^*')