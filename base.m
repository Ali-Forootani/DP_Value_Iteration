function [s] = base(m,N,n)

% We change the base of number "n" from decimal into base of "N+1"
% e.g. if "m=2", "N=2",
% if "n=1", then "s=[0 1]"
% if "n=2", then "s=[0 2]"
% if "n=3", then "s=[1 0]"

r=0;
qutient=1;
x=m;
while qutient~= 0
    qutient=floor(n/(N+1));
    remainder=mod(n,(N+1));
    s(1,x-r)=remainder;
    n=qutient;
    r=r+1;
end