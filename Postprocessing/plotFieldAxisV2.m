%plot streamslines and intensity of velocities

function plotFieldAxisV2(xMesh,yMesh,X0,Y0,uxGrid,uyGrid,linewidth,linedense,fillWhite,tilt,colorMapName,colorstream,PARAM)

Uabs = sqrt(uxGrid.^2+uyGrid.^2);

[~,h1] = contourf(X0,Y0,Uabs,500);
hold on
hStream1 = streamslice(X0,Y0,uxGrid,uyGrid,linedense);
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
hStream2 = streamslice(X0,-Y0,uxGrid,-uyGrid,linedense);
set(hStream2,'LineWidth',linewidth,'Color',colorstream)
set(h2,'LineColor','none')

plotGeometryStokesV3(xMesh,yMesh,0,[],[],[],fillWhite,PARAM);
axis([min(min(X0)) max(max(X0)) -max(max(Y0)) max(max(Y0))])
camroll(tilt)
colorbar