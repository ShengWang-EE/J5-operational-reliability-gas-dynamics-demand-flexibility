function [newmpc] = updatempcFromTransient(GFUcapacity,mpc,k,nt)
% update from gas transient result of GFU capacites and for the last step 
% electricity power flow calculation

% find which scenario is k in (j)
cumnt = cumsum([0;nt]);
j = max(find(cumnt<k));
newmpc = mpc{j};

newmpc.gen(newmpc.gfuIndex,9) = GFUcapacity';
newmpc.gen(newmpc.gfuIndex,10) = 0;


end