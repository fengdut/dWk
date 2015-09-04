function [sumy]=simpintegral(y,ny,dx)

sumy=0;
for iy=1:ny
    sumy= sumy+y(iy)*3;
end

for iy=4:3:ny-1
    sumy = sumy-y(iy);
end

    sumy = sumy-2*(y(1)+y(ny));
    sumy = dx*3*sumy/8;

end