function [Stot] = RIGHTUP_vars(Nx,Ny) 
Ntot=Nx*Ny;
Stot=zeros(Ntot,1);
for x=1:Nx
    for y=1:Ny
        i=xy2i(x,y,Nx,Ny);
        i1=xy2i(x-1,y+1,Nx,Ny); 
        Stot(i)=i1;
    end
end
end