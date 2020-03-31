%plot comparison between linear and non linear simulations at different
%time step

%set(0,'defaultaxesfontsize',20,'defaultaxeslinewidth',.7,'defaultlinelinewidth',1,'defaultpatchlinewidth',.7);

% end_loop = 200;
% time = 1:end_loop-1:end_loop;
% time = linspace(0,);
% %time = 5500:1:size(risa,2)-1;

elon_mine = zeros(PARAM.loop/PARAM.checkpoint,1);
myarea = zeros(PARAM.loop/PARAM.checkpoint,1);
V_in = axis_int(risa(:,1)',risb(:,1)');

%figure

for i = 1:PARAM.loop/PARAM.checkpoint
    
    disp(i)
    elon_mine(i) = max(risa(1:nbrel(i),i))-min(risa(1:nbrel(i),i));
    %elon_lailai(i) = max(y)-min(y);
    
    xcm = center_mass(risa(1:nbrel(i)+1,i)',risb(1:nbrel(i)+1,i)');
    
    figure(2)
    plot(risb(1:nbrel(i)+1,i),-risa(1:nbrel(i)+1,i)+xcm,'-k',-risb(1:nbrel(i)+1,i),-risa(1:nbrel(i)+1,i)+xcm,'-k','LineWidth',2)
    %plot(risa(1:nbrel(i)+1,i),risb(1:nbrel(i)+1,i),'-k',risa(1:nbrel(i)+1,i),-risb(1:nbrel(i)+1,i),'-k','LineWidth',2)
    %axis([-4.5 4.5 -1.5 1.5])
    axis([-3 3 -3 3])
    axis equal
    %title(strcat('Comparison between linear and non linear simulation for t=',num2str(time(i)),...
         %' \lambda=',num2str(visc),' Ca=',num2str(Ca)))
    %title(strcat('t=',num2str(i*PARAM.deltaT*PARAM.checkpoint)))
    title(strcat('t=',num2str(i*PARAM.deltaT*PARAM.checkpoint)))
    hold on
    %plot(xcm,0,'or')
    %plot(-x,-y,'--r',x,-y,'--r','LineWidth',2)
    hold off
    drawnow
    
    %myarea(i) = Area(loop);
    
end