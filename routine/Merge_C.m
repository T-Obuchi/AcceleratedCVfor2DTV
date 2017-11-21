function [A_mgd]=Merge_C(x,A,Active,Locked_all)
A_mgd=A;    
Ncluster=size(Locked_all,1);
for i=1:Ncluster
    Cluster=find(Locked_all(i,:));
    [var,ind]=min(abs(x(Cluster))); % Remove Representative from Cluster
    Rep=Cluster(ind);
    A_mgd(:,Rep)=sum(A(:,Cluster),2);
end
A_mgd=A_mgd(:,Active);
end