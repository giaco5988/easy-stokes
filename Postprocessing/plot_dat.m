%plot graph from initial condition.dat

close all
clear all

bef_aft = input('Before or after: ');
Ca = 6;
visc = 0.5;
Gt = input('Gt= ');
delta = input('Amplitude of the deformation= ');
color = input('Color of the graph: ');

if bef_aft==1
    cd area_norm_before
elseif bef_aft==2
    cd area_norn_max
elseif bef_aft==3
    cd area_norm_before
end

if Gt==0.1
    cd gain_10
elseif Gt==0.5
    cd gain_50
elseif Gt==0.25
    cd gain_25
end

cd(strcat('lambda=',num2str(visc),'_Ca=',num2str(Ca)))

filename = strcat('crd_rzlam',num2str(visc),'_ca',num2str(Ca),'_mA1000_delta',num2str(delta),'_time0.dat');
A = importdata(filename);
    
x = A(:,1);
y = A(:,2);

if color==1
    plot(x,y,'k',-x,y,'k','Linewidth',2)
    axis equal
elseif color==2
    plot(x,y,'r',-x,y,'r','Linewidth',2)
    axis equal
end

cd ..
cd ..
cd ..

if bef_aft==1
    saveas(gcf,strcat('Before_Gt=',num2str(Gt),'_delta=',num2str(delta),'.fig'))
elseif bef_aft==2
    saveas(gcf,stract('MAX_delta=',num2str(delta),'.fig'))
elseif bef_aft==3
    saveas(gcf,strcat('After_Gt=',num2str(Gt),'_delta=',num2str(delta),'.fig'))
end
    
    
    
    
    
    