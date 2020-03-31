% Z = peaks; surf(Z)
% axis tight
% set(gca,'nextplot','replacechildren');
% for j = 1:20
%     surf(sin(2*pi*j/20)*Z,Z)
%     F(j) = getframe;
% end
% movie(F) % Play the movie twenty times

% mov = avifile('output.avi');
% estensione = 'b.JPG';
% for i=0:239
% p = double(i);
% p = int2str(p);
% nomefile = strcat(p, estensione);
% immagine = imread(nomefile);
% f = im2frame(immagine);
% mov = addframe(mov, f);
% end;
% mov = close(mov)
fig1=figure(1);
winsize = get(fig1,'Position');
winsize(1:2)=[0 0];
numframes=800;

%movie(fig1,A,10,10,winsize)
A=moviein(numframes,fig1,winsize);
set(fig1,'NextPlot','replacechildren')
 for i=1:numframes 
     
    plot(risb(1:nbrel(i)+1,i),-risa(1:nbrel(i)+1,i),-risb(1:nbrel(i)+1,i),-risa(1:nbrel(i)+1,i),'b','LineWidth',2)
    %plot(risb(1:nbrel(i)+1,i+1),-risa(1:nbrel(i)+1,i+1),'-',-risb(1:nbrel(i)+1,i+1),-risa(1:nbrel(i)+1,i+1),'b-','LineWidth',2)
    axis equal
    ylabel('axial coordinate')
    xlabel('radial coordiante')
    title(strcat('evolving droplet t=',num2str(deltaT*i),' s'))
    xlim([-4.5 4.5])
   
   
   % add axis label, legends, titles, etc. in here 
   A(:,i)=getframe(fig1,winsize); 
 end
 
 movie(fig1,A,10,20,winsize)
 
 
 
 
 
 