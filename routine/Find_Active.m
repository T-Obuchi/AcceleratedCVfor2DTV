function [Active,Inactive]=Find_Active(x,Locked_all,Nx,Ny)
[Ncluster,Ntot]=size(Locked_all);
Belonging=[];
for i=1:Ncluster
    Cluster=find(Locked_all(i,:));
    [var,ind]=min(abs(x(Cluster))); % Remove Representative from Cluster
    Belonging=union(Belonging,[Cluster(1:ind-1),Cluster(ind+1:end)]);
end
Killed=find(abs(x)<10^(-12));       % Killed variables 
Inactive=union(Belonging,Killed);
All=[1:Ntot]';
Active=setdiff(All,Inactive);
end