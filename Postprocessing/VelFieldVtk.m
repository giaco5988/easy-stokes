%make movie from velocity field
close all

addpath('/Users/Giacomo/Documents/MATLAB/Postprocessing/writeVTK')

%movie id
ID = 2;

% Preallocate movie structure
frame = 10;
F(1:frame) = struct('cdata',[],'colormap',[]);

for i = 1:frame
    
    disp(i)
    
    [U,V,X,Y,omega] = velocity_field_buoyancy(risa,risb,risy,visc,Ca,-1.5,1.5,1.5,10,10,i+2,deltaT,checkpoint,D);
    
    X = X(:);   Y = Y(:);
    p = [X Y];
    t = delaunayn(p);
    %writeVTK(['U' num2str(i)],t,p,U);
    %writeVTK(['V' num2str(i)],t,p,V);
    writeVTK(['vorticity' num2str(i)],t,p,omega);
    F(i) = getframe; 
    close
    
end

cd ~/Desktop/
%Create AVI file.
movie2avi(F, ['VelField' num2str(ID) '.avi'], 'compression', 'None');
cd ~/Documents/MATLAB/droplet_simulations/results/rising_droplet/lambda=0.5_Ca=6/flow_field/movies/