% compute area-variation norm taking into account the eventually different
% volume (from 4/3*pi)

%clear variables
%close all

%set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',2,'defaultpatchlinewidth',.7);

path = '~/Documents/MATLAB/droplet_simulations/results/';
%load([path 'G=-1_q=100_visc=1_Dt=0.01_loop=1000_DELTA=0.2_Ca=1_RK2.mat'])

q = PARAM.q;

a = risa(1:q+1,1)';
b = risb(1:q+1,1)';

computeNow = 1;
range = 100;
step = 1;

Volume = axis_int_gauss_vect(a,b);
R = nthroot(Volume/4/pi*3,3);

%unperturbaed area
A0 = 4*pi*R^2;

if computeNow==1

    Area = zeros(range/step,1);
    V = zeros(range/step,1);

    %compute everything from geometry data
    for i=1:range/step

        display([num2str(i) ' of ' num2str(range/step)])

        aNow = risa(1:q+1,(i-1)*step+1);
        bNow = risb(1:q+1,(i-1)*step+1);

        V(i) = axis_int_gauss_vect(aNow',bNow');
        Area(i) = surf_gauss_vect(aNow',bNow');

    end

elseif computeNow==0
    
    Area = Area(1:step:range);
    V = V(1:step:range);
    V(1) = V(1)/V(1);

    
end

Area_norm = sqrt((Area-A0)/(Area(1)-A0));
%Area_norm = sqrt((Area-4*pi)/(Area(1)-4*pi));

% figure
% hold on
% if computeNow==1
% plot((0:step:range-1)*PARAM.deltaT*PARAM.checkpoint,Area)
% elseif computeNow==0
% plot((0:step:range-1)*PARAM.deltaT*PARAM.checkpoint,Area,'r')
% end
% grid on
% xlabel('t')
% ylabel('||A||')

figure
%subplot(2,1,1)
hold on
if computeNow==1
plot((0:step:range-1)*PARAM.deltaT*PARAM.checkpoint,Area-A0)
elseif computeNow==0
plot((0:step:range-1)*PARAM.deltaT*PARAM.checkpoint,Area-A0,'r')
end
grid on
xlabel('t')
ylabel('\Delta A')
title(['q=' num2str(q) ' RK=' num2str(PARAM.RK)])

figure
%subplot(2,1,1)
hold on
if computeNow==1
plot((0:step:range-1)*PARAM.deltaT*PARAM.checkpoint,Area_norm)
elseif computeNow==0
plot((0:step:range-1)*PARAM.deltaT*PARAM.checkpoint,Area_norm,'r')
end
grid on
xlabel('t')
ylabel('||A||')
title(['q=' num2str(q) ' RK=' num2str(PARAM.RK)])

figure
%subplot(2,1,2)
hold on
%try
if computeNow==1
plot((0:step:range-1)*PARAM.deltaT*PARAM.checkpoint,(V-V(1))/V(1))
elseif computeNow==0
plot((0:step:range-1)*PARAM.deltaT*PARAM.checkpoint,(V-V(1))/V(1),'r')
end
grid on
xlabel('t')
ylabel('errV')
title(['q=' num2str(q) ' RK=' num2str(PARAM.RK)])