function [a] = resultProcessing(stateVar,mpc,NK,NS,nx)
%% data dimensions
nb   = size(mpc.bus, 1);    %% number of buses
nl   = size(mpc.branch, 1); %% number of branches
ngen   = size(mpc.gen, 1);    %% number of dispatchable injections

%add GFU,LCe are included in the mpc
nGb  = size(mpc.Gbus,1); % number of gas bus
nGl  = size(mpc.Gline,1); % number of gas line
nGs  = size(mpc.Gsou,1); % number of gas source
nLCg = size(find(mpc.Gbus(:,3)~=0),1);
nEH = size(mpc.EHlocation,1);
%%
[PrsAll,GfAll,PGsAll,LCgAll,VaAll,PgAll,eiAll,giAll, ...
    eeeAll,ee3All,e13All,e1hAll,gg1All,gg2All,h3hAll, ...
    h1hAll,h14All,h24All,h2hAll,c3cAll,c4cAll,lcAll] = deal(stateVar{:});
for s = 1:NS+1
    if s == NS + 1 % 第一调度时段
        Prs1st = mat2cell(reshape(PrsAll(s,:,:),NK,sum(nx+1)),NK,nx+1);
        Gf1st = mat2cell(reshape(GfAll(s,:,:),NK,sum(nx+1)),NK,nx+1);
        [PGs1st,LCg1st,Va1st,Pg1st,ei1st,gi1st] = deal(...
            reshape(PGsAll(s,:,:),NK,nGs),reshape(LCgAll(s,:,:),NK,nLCg),...
            reshape(VaAll(s,:,:),NK,nb),reshape(PgAll(s,:,:),NK,ngen),...
            reshape(eiAll(s,:,:),NK,nEH),reshape(giAll(s,:,:),NK,nEH));
        
        
    end
        
        

end
end