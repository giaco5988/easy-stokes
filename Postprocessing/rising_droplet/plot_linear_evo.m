%plot linear evo

close all
clear variables

Ca = 6;
visc = 0.5;

set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',3,'defaultpatchlinewidth',.7,'DefaultAxesFontName','Times New Roman');
addpath(['~/Documents/MATLAB/droplet_simulations/drop_buoyancy/area_norm/lambda=' num2str(visc) '_Ca=' num2str(Ca)]);
%addpath(['~/Documents/MATLAB/droplet_simulations/results/rising_droplet/lambda=' num2str(visc) '_Ca=' num2str(Ca)]);
addpath('/Users/Giacomo/Documents/MATLAB/stability/droplet_transient/input_area_CORRECTED')

dest = '~/Documents/research_notes/APS/movie/linearEVO/';
save = 1;
AllTogether = 1;
plotEVO = 0;

%compute optimal initial shape
AAA = my_kojima(visc,Ca);     %kojima matrix
load(['finfo_Ca' num2str(Ca) '_lam' num2str(visc) '_mA1000.mat'])
F = chol(Mwgt2); %weigth matrix
%I HAVE TO "CLEAN" fopt FROM delta!!!!!!!!!!
% [~,~,V] = svd(F*expm(0.2*AAA)/F);
% fopt1 = V(:,1);
% [~,~,V] = svd(F*expm(1.05*AAA)/F);
% fopt2 = V(:,1);
% [~,~,V] = svd(F*expm(2.95*AAA)/F);
% fopt3 = V(:,1);
% [~,~,V] = svd(F*expm(5.45*AAA)/F);
% fopt4 = V(:,1);

T = time;
[~,j1] = min(abs(0.2-T));
fopt1 = fopt(:,j1);
[~,j2] = min(abs(1.05-T));
fopt2 = fopt(:,j2);
[~,j3] = min(abs(2.95-T));
fopt3 = fopt(:,j3);
[~,j4] = min(abs(5.45-T));
fopt4 = fopt(:,j4);

time = 0:0.1:10;
amp = [-0.1 -0.5 -1 0.1];
DELTA = [0.0035 0.0041 0.0056 0.0081];

%compute transient growth
% for i = 1:numel(time)
%     
%     disp(i)
%     
%     Gt(i) = norm(F*expm(5.45*AAA)/F);
%     
% end

%plot transient growth
% figure(5)
% plot(0:0.05:9.95,Gt(:,2),'k')
% grid on
% xlabel('T')
% ylabel('G_{\Delta S}')
% drawnow

DS1 = zeros(1,numel(time));
DS2 = zeros(1,numel(time));
DS3 = zeros(1,numel(time));
DS4 = zeros(1,numel(time));
%plot interface shape
for i = 1:numel(time)
        
        disp(i)
        t = time(i);
        
        if plotEVO==1
        
        %compute norm evoultion with matrix exponential
        f1 = expm(AAA*t)*fopt1;
        DS1(i) = norm(F*f1)/norm(F*fopt1);
        f2 = expm(AAA*t)*fopt2;
        DS2(i) = norm(F*f2)/norm(F*fopt2);
        f3 = expm(AAA*t)*fopt3;
        DS3(i) = norm(F*f3)/norm(F*fopt3);
        f4 = expm(AAA*t)*fopt4;
        DS4(i) = norm(F*f4)/norm(F*fopt4);
        
        %plot linear evo
        figure(5)
        plot(0:0.05:9.95,Gt(:,2),'k')
        grid on
        xlabel('T')
        ylabel('G_{\Delta S}')
        hold on
        plot(T(j1),Gt(j1,2),'om')
        plot(T(j2),Gt(j2,2),'og')
        plot(T(j3),Gt(j3,2),'ob')
        plot(T(j4),Gt(j4,2),'or')
        if i>1
        plot(time(1:i),DS1(1:i),'m')
        plot(time(i),DS1(i),'*m')
        plot(time(1:i),DS2(1:i),'g')
        plot(time(i),DS2(i),'*g')
        plot(time(1:i),DS3(1:i),'b')
        plot(time(i),DS3(i),'*b')
        plot(time(1:i),DS4(1:i),'r')
        plot(time(i),DS4(i),'*r')
        end
        hold off
        axis([0 10 0 1.2])
        %axis off
        drawnow
        
        if save==1
        name = 'LinEvo';
        print('-dpng','-loose','-r100',[dest name sprintf('%03d',i) '.png'])
        end
        
        end
        
        if AllTogether==0
        
        k = 1;
        delta = DELTA(k);
        if i==1 delta=-DELTA(k); end
        filename = strcat('crd_rzlam',num2str(visc),'_ca',num2str(Ca),'_mA1000_delta',num2str(delta),...
                '_amp',num2str(amp(k)),'_time',num2str(t),'.dat');
        A = importdata(filename);

        x = A(:,1);
        y = A(:,2);

        
%         plot(x+shift,-y+vert,this_color,'Linewidth',2)
%         plot(-x+shift,-y+vert,this_color,'Linewidth',2)
        figure(k)
        plot(x,-y,'m','Linewidth',3)
        hold on
        plot(-x,-y,'m','Linewidth',3)
        axis equal
        hold off
        axis([-1.5 1.5 -1.5 1.5])
        axis off
        drawnow
        
        if save==1
        name = 'drop1st';
        print('-dpng','-loose','-r100',[dest name sprintf('%03d',i) '.png'])
        end
        
        k = 2;
        delta = DELTA(k);
        if i==1 delta=-DELTA(k); end
        filename = strcat('crd_rzlam',num2str(visc),'_ca',num2str(Ca),'_mA1000_delta',num2str(delta),...
                '_amp',num2str(amp(k)),'_time',num2str(t),'.dat');
        A = importdata(filename);

        x = A(:,1);
        y = A(:,2);

        
%         plot(x+shift,-y+vert,this_color,'Linewidth',2)
%         plot(-x+shift,-y+vert,this_color,'Linewidth',2)
        figure(k)
        plot(x,-y,'g','Linewidth',3)
        hold on
        plot(-x,-y,'g','Linewidth',3)
        axis equal
        hold off
        axis([-1.5 1.5 -1.5 1.5])
        axis off
        drawnow
        
        if save==1
        name = 'drop2nd';
        print('-dpng','-loose','-r100',[dest name sprintf('%03d',i) '.png'])
        end
        
        k = 3;
        delta = DELTA(k);
        %if i==1 delta=-DELTA(k); end
        filename = strcat('crd_rzlam',num2str(visc),'_ca',num2str(Ca),'_mA1000_delta',num2str(delta),...
                '_amp',num2str(amp(k)),'_time',num2str(t),'.dat');
        A = importdata(filename);

        x = A(:,1);
        y = A(:,2);

        
%         plot(x+shift,-y+vert,this_color,'Linewidth',2)
%         plot(-x+shift,-y+vert,this_color,'Linewidth',2)
        figure(k)        
        plot(x,-y,'b')
        hold on
        plot(-x,-y,'b')
        axis equal
        hold off
        axis([-1.5 1.5 -1.5 1.5])
        axis off
             
        if save==1
        name = 'drop3rd';
        print('-dpng','-loose','-r100',[dest name sprintf('%03d',i) '.png'])
        end
        
        k = 4;
        delta = DELTA(k);
        if i==1 delta=-DELTA(k); end
        filename = strcat('crd_rzlam',num2str(visc),'_ca',num2str(Ca),'_mA1000_delta',num2str(delta),...
                '_amp',num2str(amp(k)),'_time',num2str(t),'.dat');
        A = importdata(filename);

        x = A(:,1);
        y = A(:,2);

        
%         plot(x+shift,-y+vert,this_color,'Linewidth',2)
%         plot(-x+shift,-y+vert,this_color,'Linewidth',2)
        figure(k)
        plot(x,-y,'r','Linewidth',3)
        hold on
        plot(-x,-y,'r','Linewidth',3)
        axis equal
        hold off
        axis off
        axis([-1.5 1.5 -1.5 1.5])
        
        if save==1
        name = 'drop4th';
        print('-dpng','-loose','-r100',[dest name sprintf('%03d',i) '.png'])
        end
        
        elseif AllTogether==1
            
            dist = 3;
            
         k = 1;
        delta = DELTA(k);
        if i==1 delta=-DELTA(k); end
        filename = strcat('crd_rzlam',num2str(visc),'_ca',num2str(Ca),'_mA1000_delta',num2str(delta),...
                '_amp',num2str(amp(k)),'_time',num2str(t),'.dat');
        A = importdata(filename);

        x = A(:,1);
        y = A(:,2);

        figure(k)
        plot(x,-y+dist,'m')
        hold on
        plot(-x,-y+dist,'m')
        axis equal
        axis off
        
        k = 2;
        delta = DELTA(k);
        if i==1 delta=-DELTA(k); end
        filename = strcat('crd_rzlam',num2str(visc),'_ca',num2str(Ca),'_mA1000_delta',num2str(delta),...
                '_amp',num2str(amp(k)),'_time',num2str(t),'.dat');
        A = importdata(filename);

        x = A(:,1);
        y = A(:,2);

        plot(x+dist,-y+dist,'g')
        plot(-x+dist,-y+dist,'g')
        axis equal
        axis off
        drawnow
        
        k = 3;
        delta = DELTA(k);
        %if i==1 delta=-DELTA(k); end
        filename = strcat('crd_rzlam',num2str(visc),'_ca',num2str(Ca),'_mA1000_delta',num2str(delta),...
                '_amp',num2str(amp(k)),'_time',num2str(t),'.dat');
        A = importdata(filename);

        x = A(:,1);
        y = A(:,2);
 
        plot(x,-y,'b')
        plot(-x,-y,'b')
        axis equal
        axis off
        
        k = 4;
        delta = DELTA(k);
        if i==1 delta=-DELTA(k); end
        filename = strcat('crd_rzlam',num2str(visc),'_ca',num2str(Ca),'_mA1000_delta',num2str(delta),...
                '_amp',num2str(amp(k)),'_time',num2str(t),'.dat');
        A = importdata(filename);

        x = A(:,1);
        y = A(:,2);

        plot(x+dist,-y,'r')
        plot(-x+dist,-y,'r')
        axis equal
        hold off
        axis off
        axis([-1.5 4.5 -1.5 4.5])
        title(['t=' num2str(t)], 'Interpreter', 'latex');
        drawnow
        
            
            if save==1
                name = 'AllTogether';
                print('-dpng','-loose','-r100',[dest name sprintf('%03d',i) '.png'])
            end
            
        end

end