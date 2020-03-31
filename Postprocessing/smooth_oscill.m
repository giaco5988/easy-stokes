%smooth oscillation for velocity field visualization

function U = smooth_oscill(Uin)

    oscill = 1;
    i = 1;

    %smooth until the criteria is respected everywhere
    %while (sum(sum(oscill)) > 0)
        
    disp(['loop' num2str(i)])

    %sum of velocities around one point
    %point in the middle
    sumCent = (Uin(1:end-2,2:end-1) + Uin(3:end,2:end-1) + Uin(2:end-1,1:end-2) + Uin(2:end-1,3:end))/4;
    %for sides
    sumUp = (Uin(1,1:end-2) + Uin(1,3:end) + Uin(2,2:end-1))/3;
    sumDown = (Uin(end,1:end-2) + Uin(end,3:end) + Uin(end-1,2:end-1))/3;
    sumLeft = (Uin(1:end-2,1) + Uin(3:end,1) + Uin(2:end-1,2))/3;
    sumRight = (Uin(1:end-2,end) + Uin(3:end,end) + Uin(2:end-1,end-1))/3;
    %for corners
    UpLeft = (Uin(1,2) + Uin(2,1))/2;
    UpRight = (Uin(1,end-1) + Uin(2,end))/2;
    DownLeft = (Uin(end,2) + Uin(end-1,1))/2;
    DownRight = (Uin(end,end-1) + Uin(end-1,end))/2;
    
    %compose sum
    sumTot = [UpLeft sumUp UpRight; sumLeft sumCent sumRight; DownLeft sumDown DownRight];
    
    %define lower and upper limit
    lower = 0.9*abs(sumTot);     upper = 1.1*abs(sumTot);
    
    %when is oscillating?
    oscill = (abs(Uin)>lower) + (abs(Uin)<upper);
    oscill = oscill==2; %0 when is oscillating
    
    Uin = Uin.*oscill;
    oscill = oscill==0; %1 when is oscillating
    U = Uin + oscill.*sumTot;
    
    Uin = U;
    
    i = i+1;
    
    %end

end