%plot streamslines and intensity of velocities

function plotFieldAxisNoMap(xMesh,yMesh,X0,Y0,uxGrid,uyGrid,linewidth,linedense,fillWhite,tilt,colorstream,PARAM)

hStream1 = streamslice(X0,Y0,uxGrid,uyGrid,linedense);
hold on
xlabel('$z$','interpreter','latex')
ylabel('$r$','interpreter','latex')
set(hStream1,'LineWidth',linewidth,'Color',colorstream)
drawnow
               
hStream2 = streamslice(X0,-Y0,uxGrid,-uyGrid,linedense);
set(hStream2,'LineWidth',linewidth,'Color',colorstream)

plotGeometryStokesV3(xMesh,yMesh,0,[],[],[],fillWhite,PARAM);
axis([min(min(X0)) max(max(X0)) -max(max(Y0)) max(max(Y0))])
camroll(tilt)
colorbar