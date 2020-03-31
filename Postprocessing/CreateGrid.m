%create grid for visualization

function [X,Y] = CreateGrid(x1,x2,y1,y2,MeshX,MeshY)

    xLenght = x2-x1;
    Nx = xLenght*MeshX;
    
    yLenght = y2-y1;
    Ny = yLenght*MeshY;

    x = linspace(x1,x2,Nx);
    y = linspace(y1,y2,Ny);
    
    [X,Y] = meshgrid(x,y);

end