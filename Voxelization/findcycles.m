function cycles=findcycles(G)
numNodes = size(G,1); 
cycles=cell(0,1);
for n = 1:numNodes
   [D,P]=graphtraverse(G,n);
   for d = D
       if G(d,n)
           cycle_temp=graphpred2path(P,d);
           if length(cycle_temp)>3
              cycles{end+1}=cycle_temp;                
           end
       end
   end
   G(n,:)=0; 
end