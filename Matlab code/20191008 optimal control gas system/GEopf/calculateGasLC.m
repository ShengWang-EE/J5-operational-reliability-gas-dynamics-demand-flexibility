function [gasLC,GFUcapacity] = calculateGasLC(transientResult,steadyResult,newmpc,nt)
%% initialization
NK = sum(nt);
% gasLoadBus = newmpc.Gbus(find(newmpc.Gbus(:,2)==2),1);
gasLoadBus = newmpc.Gbus(find(newmpc.Gbus(:,3)>0),1);
nGld = size(gasLoadBus,1);
nGb = size(newmpc.Gbus,1);
nGl = size(newmpc.Gline,1);
nGs = size(newmpc.Gsou,1);
gasLC = zeros(NK,nGld);

ngfu = size(newmpc.gfuIndex,1);
gfuUnitIndex = newmpc.gfuIndex;
GFUcapacity = zeros(NK,ngfu);
GFUelectricityBus = newmpc.gen(gfuUnitIndex,1);
GFUgasBus = zeros(size(GFUelectricityBus));
for i = 1:size(GFUelectricityBus,1)   
    GFUgasBus(i) = newmpc.GEcon(find(newmpc.GEcon(:,2)==GFUelectricityBus(i)),1); %GFU gas bus,sorted by gen order
end

gfuEfficiency = 200;


cumnt = cumsum([0;nt]);

%%
for k = 1:NK
    j = max(find(cumnt<k));
    GFUgasRequirement = steadyResult{j}.gen(newmpc.gfuIndex,2) / gfuEfficiency;
    Qinjection = zeros(nGb,1);
    % Qinjection by Q from gas pipeline
    for i = 1:nGl
        Qinjection(newmpc.Gline(i,1)) = Qinjection(newmpc.Gline(i,1)) - transientResult.Q{i}(k,1); % out-flow is assumed negtive
        Qinjection(newmpc.Gline(i,2)) = Qinjection(newmpc.Gline(i,2)) + transientResult.Q{i}(k,end);
    end
    % Qinjection by Gsou
    for i = 1:nGs
        Qinjection(newmpc.Gsou(i,1)) = Qinjection(newmpc.Gsou(i,1)) + steadyResult{j}.Gsou(i,5);
    end
    
%     % gas load of GFU
%     for i = 1:size(mpc.gfuIndex)
%         electricityBus = mpc.gen(gfuIndex(i),1);
%         gasBus = mpc.Gecon(find(mpc.GEcon(:,2)==electricityBus),1);
%         Qinjection(gasBus) = Qinjection(gasBus) + steadyResult{j}1.gen(mpc.gfuIndex(i),2) / gfuEfficiency;
%     end
%     % gas load of pure gas load
%     for i = 1:nGb
%         Qinjection(i) = Qinjection(i) + mpc.Gbus(i,3);
%     end

   % then Qinjection is the total gas injection. distribute the gas injection
   % between pure gas load and GFU
   for i = 1:nGb
       if ismember(i,gasLoadBus) && ismember(i,GFUgasBus)
           gfuIndexIndex = find(GFUgasBus==i); % the index of gfuIndex,may be a matrix
           gasLoadBusIndex = find(gasLoadBus==i);
           if Qinjection(i) <= steadyResult{j}.Gbus(i,3)
                GFUcapacity(k,gfuIndexIndex) = 0;
                gasLC(k,gasLoadBusIndex) = steadyResult{j}.Gbus(i,3) - Qinjection(i);
           end
           if Qinjection(i) > steadyResult{j}.Gbus(i,3) && Qinjection(i) < steadyResult{j}.Gbus(i,3)+sum(GFUgasRequirement(gfuIndexIndex))
               gasLC(k,gasLoadBusIndex) = 0;
               GFUcapacity(k,gfuIndexIndex) = (Qinjection(i) - steadyResult{j}.Gbus(i,3))/sum(GFUgasRequirement(gfuIndexIndex)) * GFUgasRequirement(gfuIndexIndex) * gfuEfficiency;
           end
           if Qinjection(i)>=steadyResult{j}.Gbus(i,3)+sum(GFUgasRequirement(gfuIndexIndex))
               gasLC(k,gasLoadBusIndex) = 0;
               GFUcapacity(k,gfuIndexIndex) = steadyResult{j}.gen(gfuUnitIndex(gfuIndexIndex),2);
           end
       end
       %
       if ismember(i,gasLoadBus) && ~ismember(i,GFUgasBus)
           gasLoadBusIndex = find(gasLoadBus==i);
           if steadyResult{j}.Gbus(i,3) - Qinjection(i) > 0
              gasLC(k,gasLoadBusIndex) = steadyResult{j}.Gbus(i,3) - Qinjection(i);
           else
               gasLC(k,gasLoadBusIndex) = 0;
           end
       end
       %
       if ~ismember(i,gasLoadBus) && ismember(i,GFUgasBus)
           gfuIndexIndex = find(GFUgasBus==i);
           if Qinjection(i) < sum(GFUgasRequirement(gfuIndexIndex))
                GFUcapacity(k,gfuIndexIndex) = Qinjection(i)/sum(GFUgasRequirement(gfuIndexIndex)) * GFUgasRequirement(gfuIndexIndex) * gfuEfficiency;
           else
                GFUcapacity(k,gfuIndexIndex) = GFUgasRequirement(gfuIndexIndex) * gfuEfficiency;
           end
       end
   end
end
%% force the where GFUcapacity<0 to be GFUcapacity=0
GFUcapacity(GFUcapacity<0) = 0;

%%
for i = 1:size(gasLC,2)
    gasLC(:,i) = diminishOvershootOfaCurve(gasLC(:,i));
end
for i = 1:size(GFUcapacity,2)
    GFUcapacity(:,i) = diminishOvershootOfaCurve(GFUcapacity(:,i));
end

end