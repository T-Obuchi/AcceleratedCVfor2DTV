function [H_TV]=Hessian_TV_full_delta(x,Nx,Ny,TVterms,Active_TV,delta)

Ntot=Nx*Ny;
% Geometric info
LEFTs=LEFT_vars(Nx,Ny);
RIGHTs=RIGHT_vars(Nx,Ny);
UPs=UP_vars(Nx,Ny);
DOWNs=DOWN_vars(Nx,Ny);
RIGHTUPs=RIGHTUP_vars(Nx,Ny);
LEFTDOWNs=LEFTDOWN_vars(Nx,Ny);

% Parts of Hessian in addition to TVterms
D_row_sq=Delta_Row(x,Nx,Ny);                % Difference of rows
D_col_sq=Delta_Col(x,Nx,Ny);                % Difference of columns
D_row=reshape(D_row_sq,Ntot,1);
D_col=reshape(D_col_sq,Ntot,1);

% Direct components
TVinv3=(TVterms.^2+delta^2).^(-3/2);

% Computing Hessian 
H_TV=zeros(Ntot,Ntot);
for i=1:Ntot
    Up=UPs(i);
    Down=DOWNs(i);
    Left=LEFTs(i);
    Right=RIGHTs(i);
    Leftdown=LEFTDOWNs(i);
    Rightup=RIGHTUPs(i);

    % Diagonals
    if Down~=0 & Right~=0
       H_TV(i,i)=Active_TV(i)*( (D_col(i)-D_row(i))^2+2*delta^2 )*TVinv3(i);
    end
    if Up~=0 & Rightup~=0
        H_TV(i,i)=H_TV(i,i)+Active_TV(Up)*( D_col(Up)^2+delta^2 )*TVinv3(Up);
    end
    if Left~=0 & Leftdown~=0
        H_TV(i,i)=H_TV(i,i)+Active_TV(Left)*( D_row(Left)^2+delta^2 )*TVinv3(Left);
    end

    % Nearest neighbors' 
    if Up~=0 & Rightup~=0
         H_TV(Up,i)=Active_TV(Up)*(D_row(Up)*D_col(Up)-D_col(Up)^2-delta^2)*TVinv3(Up);
    end
    if Left~=0 & Leftdown~=0
        H_TV(Left,i)=Active_TV(Left)*(D_row(Left)*D_col(Left)-D_row(Left)^2-delta^2)*TVinv3(Left);
    end
    if Down~=0 & Right~=0
        H_TV(Down,i)=Active_TV(i)*(D_row(i)*D_col(i)-D_col(i)^2-delta^2)*TVinv3(i);
    end
    if Right~=0 & Down~=0
         H_TV(Right,i)=Active_TV(i)*(D_row(i)*D_col(i)-D_row(i)^2-delta^2)*TVinv3(i);
    end

    % Next Nearest neighbors' 
    if Up~=0 & Rightup~=0 
        H_TV(Rightup,i)=Active_TV(Up)*( -D_row(Up)*D_col(Up) )*TVinv3(Up);
    end
    if Left~=0 & Leftdown~=0 
        H_TV(Leftdown,i)=Active_TV(Left)*( -D_row(Left)*D_col(Left) )*TVinv3(Left);
    end
end

end