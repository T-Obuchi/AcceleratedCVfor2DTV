function [Locked_all]=Find_Locked_all(Locked_cddt,Nx,Ny)

% Prepare variables
Ntot=Nx*Ny;
Locked_cddt2=Locked_cddt;
All=setdiff(Locked_cddt2,[0]);
Locked_all=[];                       % Memorizing all rocked cluster

while isempty(Locked_cddt2)==0
    Checked=[0];
    Pool=setdiff(Locked_cddt2(1,:),[0]); % Set the first row's candidates into the pool 
    Locked_cddt2(1,:)=[];                % Remove the first row

    % Search all the vars.. link  to vars. in Pool
    while isempty(Pool)==0
        Focused=Pool(1); Pool(1)=[];                                  % Take Focused var. from Pool 
        Checked=union(Checked,Focused);                               % Check the Focused 
        Connected=find(Locked_cddt2==Focused);                        % Find the connected vars. to Focused 
        [Sx_cddt,Sy_cddt]=size(Locked_cddt2);
        [Rows,Cols]=i2xy(Connected,Sx_cddt);                          % Column indices of Connected    
        Pool=union(Pool,sort(setdiff(Locked_cddt2(Rows,:),Checked))); % Push all vars. in the connected columns into Pool except for Checked
        Locked_cddt2(Rows,:)=[];                                      % Connected columns are removed from candidates
    end    
    Coded=zeros(1,Ntot);   
    Coded(setdiff(Checked,[0]))=1;        % Coding the indice 
    Locked_all=[Locked_all;Coded];        % Memorize all the rocked cluster 
end

end

