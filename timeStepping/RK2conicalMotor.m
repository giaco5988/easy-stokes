%Runge-Kutta order 2

function [T,X,Y,U,V] = RK2conicalMotor(f,TTT,initial,DT,initialDT,event,PARAM)

%take variables step
if initialDT==DT
    ttt = TTT(1):DT:TTT(end);
else
    ttt = [TTT(1):initialDT:TTT(20)-initialDT TTT(20):DT:TTT(end)];
end
manyDT = diff(ttt);
loop = numel(ttt);

%check dt and saving time
if sum(sum(((repmat(ttt,numel(TTT),1)>repmat(TTT',1,numel(ttt))-1e-6) + (repmat(ttt,numel(TTT),1)<repmat(TTT',1,numel(ttt))+1e-6))==2))<numel(TTT)
    warning('Saving frequency is not a correct multiple of the time step, saving steps can be lost')
end

%initialize
X = cell(numel(TTT),1);
U = cell(numel(TTT),1);
Y = cell(numel(TTT),1);
V = cell(numel(TTT),1);
T = zeros(numel(TTT),1);

%initial condition
x2 = initial{1};
y2 = initial{2};
tParametric = initial{3};

count = 1;
for i = 1:loop
    
    %current time
    t =ttt(i);
    
    %current initial condition
    x0 = x2;

    %compute solution at t0
    u0 = f(t,x0);
    
    if i~=loop
        dt = manyDT(i);

        %advance half step
        x1 = advanceBlocks(x0,dt/2,u0);
        
        %compute solution at half step
        u1 = f(t+dt/2,x1);

        %advance x0 with u1
        x2 = advanceBlocks(x1,dt,u1);
        
    end
    
    %save data
    if sum((t<TTT+dt/10) & (t>TTT-dt/10))==1
        
        %staus function
        status = event(t,x0,'');
        if status==1
            %store results
            T(count) = t;
            Y(count,:) = x0;
            V(count,:) = u0;
            T = T(1:count);
            Y = Y(1:count,:);
            V = V(1:count,:);
            break;
        end
        
        %store results
        T(count) = t;
        Y{count} = x0;
        V{count} = u0;
        
        %increase counter
        count = count+1;
        
    end

end