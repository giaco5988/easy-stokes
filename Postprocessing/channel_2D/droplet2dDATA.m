%get data from 2D droplet in a channel

clear variables
close all

%usefull path
here = pwd;
data = '~/Documents/MATLAB/droplet_simulations/results/channel_2D';

%physical variables
visc = 1;
%Ca = [0.12 0.07 0.012];
Ca = 0.012;
R = 1;  L = 25; alpha = 2.1;    Q = 1;

minH = zeros(numel(visc),numel(Ca));
constH = zeros(numel(visc),numel(Ca));
velDrop = zeros(numel(visc),numel(Ca));

for lll = 1:numel(visc)
    for rrr = 1:numel(Ca)
        
        %load data
        filename = ['Channel2D_Q=' num2str(Q) '_visc=' num2str(visc(lll)) '_Ca=' num2str(Ca(rrr)) '_R=' num2str(R) '_L=' num2str(L) '_alpha=' num2str(alpha) '_RK2.mat'];
        cd(data)
        load(filename)
        cd(here)
        %give reaches convergence
        n = PARAM.n;    m = PARAM.m;    j = PARAM.j;    q = PARAM.q;
        ind = find(risa(2,:)==0,1,'first');
        if isempty(ind)
            ind = PARAM.loop/PARAM.checkpoint;
        end
        aDrop = risa(n+2*m+j+2:end,ind-1);  bDrop = risb(n+2*m+j+2:end,ind-1);
        
        %Uavg of the inlet
        Uavg = PARAM.Q/2/PARAM.R;

        %minimum film thickness
        minH(lll,rrr) = min(max(risb(:,ind-1))-bDrop)/2/PARAM.R;
        constH(lll,rrr) = min(max(risb(:,ind-1))-bDrop(q/4))/2/PARAM.R;
        velDrop(lll,rrr) = risy(2*(n+2*m+j)+1,ind-2)/Uavg;
        
    end
end

figure
subplot(3,1,1)
plot(Ca,minH,'o-')
grid on
xlabel('Ca')
ylabel('min(h)')
title(['\lambda=' num2str(visc)])

subplot(3,1,2)
plot(Ca,constH,'o-')
grid on
xlabel('Ca')
ylabel('h')
title(['\lambda=' num2str(visc)])

subplot(3,1,3)
plot(Ca,velDrop,'o-')
grid on
xlabel('Ca')
ylabel('v_{drop}')
title(['\lambda=' num2str(visc)])














