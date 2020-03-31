%compute concentration field

function [Xsing,Ysing,ux,uy,p] = computeVelPressField(X,Y,x,y,solution,velBlocks,Uback,PARAM,substitutePoint,coeffSub,noIn)

% if noIn==1
%     warning('Sometimes it can provoke inaccurate substituition of boundary points')
% end

%change variables name
PARAM.typeBC = PARAM.typeBCstokes;
PARAM.orderVariable = PARAM.orderVariableStokes;
PARAM.orderGeometry = PARAM.orderGeometryStokes;

%compute number of element (sounting constant and linear differently)
[xStore,yStore,nnn] = computeSingularityLocation(x,y,PARAM);

[line,column] = size(X);
Xsing = X(:);
Ysing = Y(:);

%compute single and double layer operator
[GXX,GXY,GYX,GYY,A11,A12,A21,A22,PX,PY,PI1,PI2] = computeKernelsOperatorsAxisVisu(Xsing,Ysing,x,y,PARAM);

%assembly single and double layer operator
[SL,DL,SLpressure,DLpressure,fBC,uBC] = assemblySingkeAndDoubleLayer(solution,GXX,GXY,GYX,GYY,A11,A12,A21,A22,PX,PY,PI1,PI2,Xsing,nnn,x,y,PARAM);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure out if I'm inside or outside a block
InOutBefore = 0;
maxElem = 0;

if noIn==1
    for i = 1:numel(PARAM.panels)

        %compute block coordinates
        [xHere,yHere,~,~,~,dlHere] = getBlockCoordinates(x,y,PARAM,i);

        %compute maximum elem
        maxElem = max(max(dlHere),maxElem);

        InOut = FigInOut(xHere,yHere,Xsing,Ysing);
        InOut = (InOut>pi)+InOutBefore;

        InOutBefore = InOut;

    end
    InOut = (InOut>0);
else
    InOut = 0;
end

%add repulsive force for droplet (or bubble)
for i = 1:numel(PARAM.panels)
      
      if (PARAM.repulsiveForces(i)==3||PARAM.repulsiveForces(i)==4) && PARAM.blockType(i)==2
      
      if numel(PARAM.panels)>2
          error('Current implementation works only for two block')
      end
      
      %range
      [~,~,~,~,panelRange] = getBlockCoordinates(x,y,PARAM,i);
      startMatrix = getSingRange(panelRange(1),nnn);
      [~,endMatrix] = getSingRange(panelRange(end),nnn);
      
      for k = 1:numel(PARAM.panels)
          
            if k~=i
      
                [dfXrep,dfYrep] = disjoiningPressureBlocks(x,y,i,k,PARAM);
                fBC(2*startMatrix:2:2*endMatrix) = fBC(2*startMatrix:2:2*endMatrix) + dfYrep;
                if PARAM.repulsiveForces(i)==4
                   fBC(2*startMatrix-1:2:2*endMatrix-1) = fBC(2*startMatrix-1:2:2*endMatrix-1) + dfXrep;
                end
            
            end
            
      end
    
      end
    
end

% find the velocity field solving Ufield = -SL*f_bc + DL*u_bc
u = - SL*fBC/8/pi + DL*uBC/8/pi;
ux = u(1:2:end-1);
uy = u(2:2:end);

% find the velocity field solving Pfield = -SL*f_bc + DL*u_bc
p = - SLpressure*fBC/8/pi + DLpressure*uBC/8/pi;

%check if the point is inside a droplet
for i = 1:numel(PARAM.panels)
    
    if PARAM.blockType(i)==2

        %compute block coordinates
        [xHere,yHere,~,~,panelRange,dlHere] = getBlockCoordinates(x,y,PARAM,i);

        %compute maximum elem
        maxElem = max(max(dlHere),maxElem);

        InDrop = FigInOut(xHere,yHere,Xsing,Ysing);
        InDrop = (InDrop>pi);

        visc = PARAM.visc(panelRange);
        if visc==0
            ux = ux.*(1-InDrop);
            uy = uy.*(1-InDrop);
            p = p.*(1-InDrop);
        else
            ux = ux.*(1-InDrop) + ux.*InDrop/visc;
            uy = uy.*(1-InDrop) + uy.*InDrop/visc;
            p = p.*(1-InDrop) + p.*InDrop/visc;
        end
        
    end

end

%add flow
if PARAM.addFlow==1
    ux = Uback + ux;
elseif PARAM.addFlow==2
    ux = Uback(1:numel(X)) + ux;
    uy = Uback(numel(X)+1:end) + uy;
end

%check if evaluation point is very close to the interface and change it if
%it is so
if substitutePoint==1
    
    [Xsing,Ysing,ux,uy,p,InOut] = closeInterfaceVisuBlocks(x,y,xStore,yStore,Xsing,Ysing,ux,uy,p,InOut,uBC,maxElem,coeffSub,velBlocks,PARAM);
    
end

%InOut = 0;
% ux = ux.*(1-InOut');
% uy = uy.*(1-InOut');
% p = p.*(1-InOut');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Xsing = reshape(Xsing,line,column);
Ysing = reshape(Ysing,line,column);
ux = reshape(ux,line,column);
uy = reshape(uy,line,column);
p = reshape(p,line,column);








