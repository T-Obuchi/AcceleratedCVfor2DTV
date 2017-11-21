function [D_col] = Delta_Col(x,Nx,Ny)
Ntot=Nx*Ny;
x_sq = reshape(x,Nx,Ny);

% x: vertical direction, y: horizontal direction
D_col=x_sq(:,2:Ny)-x_sq(:,1:Ny-1);                   % columns' difference
D_col=[D_col,zeros(Nx,1)];
end
