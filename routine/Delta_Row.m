function [D_row] = Delta_Row(x,Nx,Ny)
Ntot=Nx*Ny;
x_sq = reshape(x,Nx,Ny);

% x: vertical direction, y: horizontal direction
D_row=x_sq(2:Nx,:)-x_sq(1:Nx-1,:);                   % rows' difference
D_row=[D_row;zeros(1,Ny)];
end
