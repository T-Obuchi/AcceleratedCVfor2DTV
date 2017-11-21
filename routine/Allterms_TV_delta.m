function [TVterms]=Allterms_TV_delta(x,Nx,Ny,delta)
Ntot=Nx*Ny;
x_sq = reshape(x,Nx,Ny);

D_row=Delta_Row(x,Nx,Ny);
D_col=Delta_Col(x,Nx,Ny);
TVterms_sq=sqrt(D_col.^2+D_row.^2+delta^2); % All terms of TV in matrix form
TVterms=reshape(TVterms_sq,Ntot,1); % All terms of TV in verctor form
end
