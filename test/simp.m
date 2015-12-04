function f=simp(n,h,fi)
s=0.0;
s0=0.0;
s1=0.0;
s2=0.0;
for j=2:2:(n-1)
    s1=s1+fi(j-1);
    s0=s0+fi(j);
    s2=s2+fi(j+1);
end
s=h*(s1+4.0*s0+s2)/3.0;

if mod(n,2)==0.0
    s=s+h*(5.0*fi(n)+8.0*fi(n-1)-fi(n-2))/12.0;
end
f=s;