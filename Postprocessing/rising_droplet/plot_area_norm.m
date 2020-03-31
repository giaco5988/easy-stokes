%plot area-variation norm at certain times

function plot_area_norm(time,delta)

    Area = zeros(numel(time),1);

    for k=1:numel(time)

        %disp(k)
        filename = strcat('crd_rzlam0.5_ca6_mA1000_delta',num2str(delta),'_time',num2str(time(k)),'.dat');
        A = importdata(filename);

        x = A(:,1);
        y = A(:,2);

        %V2 = axis_int_gauss(y',x');
        Area(k) = surf_gauss(y',x');
        
    end
    
    figure
    norm = sqrt((Area-4*pi)/(Area(1)-4*pi));
    plot(time,norm,'o-')
    xlabel('time')
    ylabel('A norm')
    grid on
    
    %savefig(strcat('a_norm_D=',num2str(delta),'.fig'))
    %close
    
    figure(1)
    hold all
    plot(time,norm,'o-')
    xlabel('time')
    ylabel('A norm')
    hold off
    grid on
    
end