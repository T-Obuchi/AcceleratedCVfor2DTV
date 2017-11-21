function [i] = xy2i(x,y,Nx,Ny)
if x<1 || x>Nx || y<1  || y>Ny  
    i=0;
else
    i=x+Nx*(y-1);
end
end

