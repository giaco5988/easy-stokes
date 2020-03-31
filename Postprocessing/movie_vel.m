%make movie from velocity field
close all

%movie id
ID = 2;

% Preallocate movie structure
frame = 250;
F(1:frame) = struct('cdata',[],'colormap',[]);

for i = 1:frame
    
    disp(i)
    
    [U,V,X,Y] = velocity_field_buoyancy(risa,risb,risy,visc,Ca,-1.5,1.5,1.5,10,10,i+2,deltaT,checkpoint,D);
    
    scrsz = get(0,'ScreenSize');
    %figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2])
    
    figure(1)
    set(gcf,'Position',[1 scrsz(4)/1 scrsz(3)/1 scrsz(4)/1])
    saveas(gcf,['~/Documents/research_notes/DTU_summer_school/presentation/movie1/quiver' num2str(i) '.png'])
    close
    
    figure(2)
    set(gcf,'Position',[1 scrsz(4)/1 scrsz(3)/1 scrsz(4)/1])
    saveas(gcf,['~/Documents/research_notes/DTU_summer_school/presentation/movie1/vorticity' num2str(i) '.png'])
    close
    
    figure(3)
    set(gcf,'Position',[1 scrsz(4)/1 scrsz(3)/1 scrsz(4)/1])
    saveas(gcf,['~/Documents/research_notes/DTU_summer_school/presentation/movie1/interface' num2str(i) '.png'])
    close
    
end

%cd ~/Desktop/
%Create AVI file.
%movie2avi(F, ['VelField' num2str(ID) '.avi'], 'compression', 'None');
%cd ~/Documents/MATLAB/droplet_simulations/results/rising_droplet/lambda=0.5_Ca=6/flow_field/ellipsoidal