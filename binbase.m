function [s] = binbase(m,n)

% We change the base of number "n" from decimal into binary base.
% e.g. binbase(3,5)=[1 0 1]
% e.g. binbase(2,1)=[0 1]
% e.g. binbase(3,1)=[0 0 1]
% e.g. binbase(3,3)=[0 1 1]
% e.g. binbase(3,4)=[1 0 0]

r=0;
qutient=1;
x=m;
while qutient~= 0
    qutient=floor(n/(2));
    remainder=mod(n,(2));
    s(1,x-r)=remainder;
    n=qutient;
    r=r+1;
end