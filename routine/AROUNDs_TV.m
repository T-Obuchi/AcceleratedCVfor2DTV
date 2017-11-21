function [Stot] = AROUNDs_TV(Nx,Ny) 
% Nearest neighbors are defined in matching to TV convention

Ntot=Nx*Ny;
Stot=zeros(Ntot,2);
for x=1:Nx
    for y=1:Ny
        i=xy2i(x,y,Nx,Ny);
%        i1=xy2i(x-1,y,Nx,Ny); i2=xy2i(x,y-1,Nx,Ny);
        i3=xy2i(x+1,y,Nx,Ny); i4=xy2i(x,y+1,Nx,Ny);
        Stot(i,:)=[i3,i4];       
    end
end

end