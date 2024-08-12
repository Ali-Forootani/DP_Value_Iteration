function [transition combined_probability] = stateanalysis(s,N,m,a,lambda,mu)
%% This function takes a state "s", number of resources "N",
% number of prices "m", the index to one of the "m" prices "i=a",
% "lambda's", and "mu's".
% 
%
% Outputs: matrix "transition" and matrix "combined_probability".

% matrix "transition": each row of this matrix represents a possible transition
% from current state "s=x(k)" to the next state "x(k+1)".
% 
% matrix "combined_probability": corresponding to each row of "transition"
% matrix there is a probability of going from state "s=x(k)" to the next
% state "x(k+1)".
%===========================================

% Given a state, if an action "i=a" is taken we have to find out how many
% transitions are possible?

% We clasify our states into different types! 
% In general we have 11 type states if we consider transitions on the boundary.
% For each type we can find a closed formula regarding number of
% transitions if a specific action has been taken.
%===========================================

% "binbase": It is a written function which creates the next state "x(k+1)"
% from "s=x(k)", e.g. if "x(k)=[1 1]" then "x(k+1)" could be [1 1], [0 0], [1 0], [0 1], [2 0]
% if action "i=1" has been taken
% 
% "expectation": It is a written function which gives the probability of a
% transtion from a state "s=x(k)" to the state "x(k+1)". This function constructs
% each row of "combined_probability" matrix.

%===========================================
% auxiliary variables
n2=1;
NZ=0;
n3=0;
n5=0;
n6=0;
n4=0;
n7=0;
n8=0;
n9=0;
n10=0;
n11=0;
n12=0;
n13=0;
n14=0;
n15=1;
n16=0;
n17=1;
n18=0;
n19=1;
n20=0;
n21=1;
n22=0;
n23=1;
n24=0;
BZ=0;
DZ=0;
OR=0;
clear transition tr
%% Suppose that state "s=[1 0 0 3 4]"
% "NZ": number of zero elements in a vector "s", e.g. above NZ=2
% and there is no customer corresponding to price "c_2" and "c_3".

% "DZ": number of non-zero elements in a vector "s", e.g. for the above
% example BZ=3 and there is 1 customer in price "c_1", 3 customers in price
% "c_4", and 4 customers in price "c_5"

% "BZ": number of non-zero elements in a vector "s", and non-zero elements
% are not in the group related to action "c_i=a"

% "OR": number of zero elements in a vector "s", and zero elements are not
% in the group related to action "c_i=a"

while n2 <= m
    if s(n2)==0
        NZ=NZ+1;
        n3=n3+1;
        inNZ(1,n3)=n2;
    end
    if s(n2)~=0 
        n4=n4+1;
        DZ=DZ+1;
        inDZ(1,n4)=n2;
    end
    if s(n2)~=0 & n2~= a
        n9=n9+1;
        BZ=BZ+1;
        inBZ(1,n9)=n2;
    end
    if s(n2)==0 & n2~= a
        n24=n24+1;
        OR=OR+1;
        inOR(1,n24)=n2;
    end
    n2=n2+1;
end

% Having to know how many non-zero elements, zero elements, zero customer
% corresponding to the "c_i=a" help us to figure out for a given state how
% many transitions could be possible
% 
% when we have found the number of transitions, we have to find the
% corresponding states

%%

%In general, we have two conditions: being on the "boundry" or not!
% so at first we check "sum(s)<N" 

% "type=1": state like s=[1 2 3], N=7, m=3, and action "a=c_2"
% e.g. here we have 3*2^(2) transitions

%===============

% "type=2": state like s=[1 2 0 2], N=7, m=4, and action "a=c_2"
% e.g. here we have 3*2^(2)=12 transitions, since "x_3=0"!
%===============

% "type=3": state like s=[3 0 0 0], N=7, m=4, action "a=c_1"
% e.g. here we have 3 transitions, since "x_2=0,x_3=0,x_4=0"!

%===============
% "type=4": state like s=[0 0 0 0], N=7, m=4, action "a=c_1"
% e.g. here we have 2 transitions, since we are at the origin

%===============
% "type=5": state like s=[0 2 1 1], N=7, m=4, action "a=c_1", "x_1=0"!!!!
% e.g. here we have 2*2^(4-1)=16, since "x_1=0" and we choose "a=c_1"

% ==============
% "binbase": the aim of this function is to calculate an "increament",
% or a "decreament" in the number of resources, since we assumed only one
% increament or one decreament at each group, e.g. if we were at the state
% x=[1 1], N=3, m=2, and "a=c_1" then we have possible transitions [1 1],[2 1]
% [0 1],[2 0],[1 0],[0 0]; we are allowed 3 possibilites of "increament",
% "decreament" and "nothing" for price "c_1", and 2 possibilities of
% "decreament" and "nothing" for price "c_2".
% We do this procedure by using binary value, e.g. for above example 
% an "increament" regarding "c_1" is [1 0], a "decreament" regarding "c_2" is
% [-1 0], and "nothing" is [0 0].
% ==============





if sum(s)<N
    if NZ==0
        transition=zeros(3*2^(m-1-NZ),m);
        type=1;
        while n5 <= 2^m-1
            tr=-1*binbase(m,n5);
            n5=n5+1;
            transition(n5,1:end)=s+tr;
        end
        while n6 <= 2^(m-1)-1
            tr=-1*binbase(m-1,n6);
            n6=n6+1;
            n5=n5+1;
            transition(n5,a)=s(a)+1;
            transition(n5,inBZ)=s(1,inBZ)+tr;
        end 
    end
    if NZ~=0 & s(a)~=0 & NZ<(m-1)
            transition=zeros(3*2^(m-1-NZ),m);
            type=2;
            while n7 <= 2^(m-NZ)-1
                tr=-1*binbase(m-NZ,n7);
                n7=n7+1;
                transition(n7,inDZ)=s(1,inDZ)+tr;
            end
            while n8 <= 2^(m-NZ-1)-1
                tr=-1*binbase(m-NZ-1,n8);
                n7=n7+1;
                n8=n8+1;
                transition(n7,inBZ)=s(inBZ)+0;
                transition(n7,a)=s(a)+1;
                transition(n7,inBZ)=s(1,inBZ)+tr;
            end
    end
    if NZ~=0 & s(a)~=0 & NZ==(m-1)
              transition=zeros(3*2^(m-1-NZ),m);
              type=3;
              while n10 <= 2^(m-NZ)-1
                tr=-1*binbase(m-NZ,n10);
                n10=n10+1;
                transition(n10,inDZ)=s(1,inDZ)+tr;
              end
              while n11 <= 2^(m-NZ-1)-1
                n10=n10+1;
                transition(n10,a)=s(a)+1;
                n11=n11+1;
              end
    end
    if s==0
                type=4;
                transition=zeros(2,m);
                transition(1,1:end)=s;
                transition(2,a)=s(a)+1;
    end
    
    if s(a)==0 & BZ~=0
        transition=zeros(2*2^(BZ),m);
        type=5;
        while n12 <= 2^(BZ)-1
            tr=-1*binbase(BZ,n12);
            n12=n12+1;
            transition(n12,inDZ)=s(1,inDZ)+tr;
        end
        while n13 <= 2^(BZ)-1
            tr=-1*binbase(BZ,n13);  
            n13=n13+1;
            n12=n12+1;
            transition(n12,inDZ)=s(1,inDZ)+tr;
            transition(n12,a)=s(a)+1;
        end
    end
    type;
    combined_probability=expectation(s,N,m,a,transition,lambda,mu,type);
end

%% If "sum(s)=N" on the Boundry
% "type=6", state like s=[1 1], N=2, m=2, "a=c_2"
% here we have 5 transitions

%===========================
% "type=7", state like s=[1 2 3 0], N=6, m=4, "a=c_1" 
% e.g. here we have 11 transitions
%===========================
% "type=8", state like s=[6 0 0 0], N=6, m=4, "a=c_1"
% e.g. here we have 2 transitions, since "x_2=0,x_3=0,x_4=0"
%===========================
% "type=9", state like s=[5 0 1 0], N=6, m=4, "a=c_2"
% e.g. here we have 7 transitions
%===========================
% "type=10", state like s=[6 0 0 0], N=6, m=4, "a=c_2"
% e.g. here we have 3 transitions
%===========================
% "type=11", state like s=[0 8 0 0], N=8, m=4, "a=c_2"
% e.g. here we have 2 transitions
% 

if sum(s)==N
    if NZ==0
        transition=zeros(1,m);
        type=6;
        while n14 <= 2^m-1
           tr=-1*binbase(m,n14);
           n14=n14+1;
           transition(n14,1:end)=s+tr;
        end
         while n15 <= 2^(m-1)-1
            tr=-1*binbase(BZ,n15);
            n15=n15+1;
            if (sum(tr)+(s(a)+1))<= N
                n14=n14+1;
                transition(n14,inBZ)=s(1,inBZ)+tr;
                transition(n14,a)=s(a)+1;
            end
         end
    end
    
    if NZ~=0 & s(a)~=0 & NZ<(m-1)
        transition=zeros(1,m);
        type=7;
        while n16 <= 2^(m-NZ)-1
           tr=-1*binbase(DZ,n16);
           n16=n16+1;
           transition(n16,inDZ)=s(inDZ)+tr; 
        end
         while n17 <= 2^(m-NZ-1)-1
            tr=-1*binbase(BZ,n17);
            n17=n17+1;
            if (sum(tr)+(s(a)+1))<= N
                n16=n16+1;
                transition(n16,inBZ)=s(1,inBZ)+tr;
                transition(n16,a)=s(a)+1;
            end
         end
    end
    
    if NZ~=0 & s(a)~=0 & NZ==(m-1)
        transition=zeros(1,m);
        type=8;
        while n18 <= 2^(m-NZ)-1
           tr=-1*binbase(DZ,n18);
           n18=n18+1;
           transition(n18,inDZ)=s(inDZ)+tr; 
        end
        while n19 <= 2^(m-NZ-1)-1
            n19=n19+1;
            if (sum(tr)+(s(a)+1))<= N
                n18=n18+1;
                transition(n18,a)=s(a)+1;
            end
        end
    end

    if NZ~=0 & s(a)==0 & NZ< (m-1)
        transition=zeros(1,m);
        type=9;
        while n20 <= 2^(BZ)-1
           tr=-1*binbase(DZ,n20);
           n20=n20+1;
           transition(n20,inDZ)=s(inDZ)+tr; 
        end
        while n21 <= 2^(BZ)-1
            tr=-1*binbase(BZ,n21);  
            n21=n21+1;
            n20=n20+1;
            transition(n20,inDZ)=s(1,inDZ)+tr;
            transition(n20,a)=s(a)+1; 
        end
    end
    
    if s(a)==0 & NZ == (m-1)
        transition=zeros(2^(m-NZ)+(m-NZ),m);
        type=10;
        while n22 <= 2^(BZ)-1
            tr=-1*binbase(BZ,n22);
            n22=n22+1;
            transition(n22,inBZ)=s(inBZ)+tr;
        end
        while n23 <= 2^(BZ)-1
            n22=n22+1;
            tr=-1*binbase(BZ,n23);  
            n23=n23+1;
            transition(n22,inBZ)=s(1,inBZ)+tr;
            transition(n22,a)=s(a)+1; 
        end
    end
    
    if s(a)==N 
        type=11;
        transition=zeros(2,m);
        transition(1,1:end)=s;
        transition(2,a)=s(a)-1;
    end
    type;
        combined_probability=expectation(s,N,m,a,transition,lambda,mu,type);
end
