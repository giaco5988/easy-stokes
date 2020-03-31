%compute concentration field

function [Xsing,Ysing,PHIfield] = computeConcentrationField(X,Y,x,y,yOut,PARAM,substituiteConc,coeffSub)

%change variables name
PARAM.typeBC = PARAM.typeBClaplace;
PARAM.orderVariable = PARAM.orderVariableLaplace;
PARAM.orderGeometry = PARAM.orderGeometryLaplace;

%compute number of element (sounting constant and linear differently)
[xStore,yStore,nnn] = computeSingularityLocation(x,y,PARAM);

[line,column] = size(X);
Xsing = X(:);
Ysing = Y(:);
%preallocation
G1 = cell(numel(PARAM.n),1);
G2 = cell(numel(PARAM.n),1);
for i = 1:numel(PARAM.n)
            
            if PARAM.orderVariable(i)==0 && PARAM.orderGeometry(i)==0   % integration on constant, straight element
        
                %compute with MATLAB  
                [G1{i},G2{i}] = computeGaxisNoSing(x{i},y{i},Xsing,Ysing);
        
            elseif PARAM.orderVariable==1 && PARAM.orderGeometry==0     % integration on linear, straight element

                error('Not implemented')

            elseif PARAM.orderVariable==0 && PARAM.orderGeometry==1     % integration on constant, curved element

                error('Not implemented')

            elseif PARAM.orderVariable==1 && PARAM.orderGeometry==1     % integration on linear, curved element

                error('Not implemented')

            end
        
end

%build matrices for single layer and double layer potential
SL = zeros(numel(Xsing),sum(nnn));
DL = zeros(numel(Xsing),sum(nnn));
PHIbc = zeros(sum(nnn),1);
gradPHIbc_n = zeros(sum(nnn),1);
for i = 1:numel(PARAM.n)
      
      if i==1
          startMatrix = 1;
      else
          startMatrix = 1+sum(nnn(1:i-1));
      end
      
      %single layer
      SL(:,startMatrix:sum(nnn(1:i))) = G1{i};
      
      %double layer
      DL(:,startMatrix:sum(nnn(1:i))) = G2{i};
      
      if PARAM.typeBC(i)==1    %prescribe phi and find gradient
          
          %known BC
          PHIbc(startMatrix:sum(nnn(1:i))) = ones(PARAM.n(i),1).*PARAM.concBC{i};
          
          %BC from the result
          gradPHIbc_n(startMatrix:sum(nnn(1:i))) = yOut(startMatrix:sum(nnn(1:i)));
          
      elseif PARAM.typeBC(i)==2    %prescribe gradient and find phi
          
          %known BC
          gradPHIbc_n(startMatrix:sum(nnn(1:i))) = ones(PARAM.n(i),1).*PARAM.fluxBC{i};
          
          %BC from the result
          PHIbc(startMatrix:sum(nnn(1:i))) = yOut(startMatrix:sum(nnn(1:i)));
          
      end
      
end

% find the concentration field solving PHIfield = -SL*gradPHIbc_n + DL*PHIbc
PHIfield = -SL*gradPHIbc_n + DL*PHIbc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure out if I'm inside or outside a block
InOutBefore = 0;
maxElem = 0;

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

%check if evaluation point is very close to the interface and change it if
%it is so
if substituiteConc==1
        [Xsing,Ysing,PHIfield] = close_interface_visuLaplace(xStore,yStore,Xsing,Ysing,PHIfield,InOut,PHIbc,maxElem,coeffSub);
end

%deactivate points inside blocks
PHIfield = PHIfield.*(1-InOut);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Xsing = reshape(Xsing,line,column);
Ysing = reshape(Ysing,line,column);
PHIfield = reshape(PHIfield,line,column);








