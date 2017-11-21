function [x,y] = i2xy(i,Nx)
y=idivide(int32(i-1),Nx)+1;
x=int32(i)-Nx*(y-1);
end

