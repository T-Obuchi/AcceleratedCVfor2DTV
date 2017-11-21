function [H_mgd]=Merge_RC(x,H,Active,Locked_all)
[Ncluster,Ntot]=size(Locked_all);
H_copy=H;
for i=1:Ncluster
    Cluster=find(Locked_all(i,:));
    [var,ind]=min(abs(x(Cluster))); 
    Rep=Cluster(ind);                  % Representative
    H_copy(:,Rep)=sum(H_copy(:,Cluster),2);
    H_copy(Rep,:)=sum(H_copy(Cluster,:),1);
end
H_mgd=H_copy(Active,Active);
end