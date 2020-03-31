%plot streamslines and intensity of velocities

function plotFieldAxisStreamlines(xMesh,yMesh,X0,Y0,uxGrid,uyGrid,linewidth,fillWhite,tilt,colorMapName,colorstream,startX,startY,yStokes,PARAM)
                       
startX1 = startX{1};
startX2 = startX{2};
startX3 = startX{3};
startY1 = startY{1};
startY2 = startY{2};
startY3 = startY{3};

Uabs = sqrt(uxGrid.^2+uyGrid.^2);

[~,h1] = contourf(X0,Y0,Uabs,500);
hold on
hStream1 = streamline(X0,Y0,uxGrid,uyGrid,[startX1 startX2 startX3],[startY1 startY2 startY3]);
%hStream1 = streamslice(X0,Y0,uxGrid,uyGrid,linedense);
set(h1,'LineColor','none')
%axis equal
xlabel('$z$','interpreter','latex')
ylabel('$r$','interpreter','latex')
colormap(colorMapName)
set(hStream1,'LineWidth',linewidth,'Color',colorstream)
%caxis([0 max(max(Uabs))])
%axis off
drawnow
               
[~,h2] = contourf(X0,-Y0,Uabs,500);
hStream2 = streamline(X0,-Y0,uxGrid,-uyGrid,[startX1 startX2 startX3],-[startY1 startY2 startY3]);
%hStream2 = streamslice(X0,-Y0,uxGrid,-uyGrid,linedense);
set(hStream2,'LineWidth',linewidth,'Color',colorstream)
set(h2,'LineColor','none')

plotGeometryStokesV3(xMesh,yMesh,0,[],[],[],fillWhite,PARAM);
axis([min(min(X0)) max(max(X0)) -max(max(Y0)) max(max(Y0))])
camroll(tilt)
colorbar
axis equal

% xyCoord1 = stream2(X0,Y0,uxGrid,uyGrid,startX1,startY1);
% xyCoord2 = stream2(X0,Y0,uxGrid,uyGrid,startX2,startY2);
% xyCoord3 = stream2(X0,Y0,uxGrid,uyGrid,startX3,startY3);
% 
% xyCoord = {xyCoord1 xyCoord2 xyCoord3};
xyCoord = stream2(X0,Y0,uxGrid,uyGrid,[startX1 startX2 startX3],[startY1 startY2 startY3]);

for zzz = 1:numel(xyCoord)
                       
%                       if isempty(xyCoord{3})==0
%                       xyCoordHere = xyCoord{zzz};
%                       xyCoordHere = xyCoordHere{1};
                      xyCoordHere = xyCoord{zzz};
                      %xCoord = xyCoordHere(50:240:end,1);
                      %yCoord = xyCoordHere(50:240:end,2);
                      sizeHere = size(xyCoordHere);
                      xCoord = xyCoordHere(ceil(sizeHere(1)/3),1);
                      yCoord = xyCoordHere(ceil(sizeHere(1)/3),2);
                      
                      %compute velocity
                      [~,~,uxLine,uyLine] = computeVelPressField(xCoord,yCoord,xMesh,yMesh,yStokes,yStokes(end),0,PARAM,0,nan,0);
                      
                      %quiver(xCoord,yCoord,uxLine,uyLine,'r','LineWidth',3)
                      %scaleHere = 5;
                      %quiver(xCoord,yCoord,scaleHere*uxLine./sqrt(uxLine.^2+uyLine.^2),scaleHere*uyLine./sqrt(uxLine.^2+uyLine.^2),'r','LineWidth',3,'AutoScale','Off','MaxHeadSize',10)
                      
                      dx = uxLine./sqrt(uxLine.^2+uyLine.^2)/100;
                      dy = uyLine./sqrt(uxLine.^2+uyLine.^2)/100;
                      dlHere = sqrt(diff(xCoord).^2+diff(yCoord).^2);
                      distOrigin = sqrt(xCoord.^2+yCoord.^2);
                      for sss = 1:numel(xCoord)
                          
                          if sss==1 && distOrigin(sss)>0
                              arrow([xCoord(sss) yCoord(sss)],[xCoord(sss)+dx(sss) yCoord(sss)+dy(sss)],'EdgeColor','k','FaceColor','w')
                              arrow([xCoord(sss) -yCoord(sss)],[xCoord(sss)+dx(sss) -yCoord(sss)-dy(sss)],'EdgeColor','k','FaceColor','w')
                          end
                          
                          if sss>1
                              if dlHere(sss-1)>40 && distOrigin(sss)>0
                                arrow([xCoord(sss) yCoord(sss)],[xCoord(sss)+dx(sss) yCoord(sss)+dy(sss)],'EdgeColor','k','FaceColor','r')
                                arrow([xCoord(sss) -yCoord(sss)],[xCoord(sss)+dx(sss) -yCoord(sss)-dy(sss)],'EdgeColor','k','FaceColor','r')
                              end
                          end
                      
                      end
                      %end
                       
end