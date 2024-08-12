function [combined_probability] = expectation(s,N,m,a,transition,lambda,mu,type)

%% In this fuction we want to know the probability of going from 
%  the state "x(k)" to "x(k+1)" based on matrix "transition" of
%  specific action "a=c_i"

%  e.g. if x(k)=[1 2 1], one row of matrix "transition" represents 
%  a possible transition x(k+1)=[1 1 1] and "c_1=a";
%  a row of matrix "combined_probability"
%  probability will be "[1-lambda_1-mu_1  mu_2  1-lambda_3]" and 
% (1-lambda_1-mu_1)*(mu_2)*(1-lambda_3) will be the probability of going
% from state [1 2 1] to state [1 1 1] is action "c_1" is taken.

%  "type": {1,2,...,11} from function "state_analysis"
%  "lambda's" and "mu's" are the probabilities corresponding to each group.
%  "N": number of resources.
%  "m": number of prices.
%  "transition": matrix produced by "state_analysis" function, where each
%  row represents a possible state "x(k+1)".

% We know that the sumamtion of all the probabilities has to be 1, so we check
% it at the end of this function the name of variable is called "summation".



%%
NT=size(transition);
n1=1;
n2=1;
n3=0;
n4=0;
n5=0;
n6=1;
n7=0;
n8=0;
n9=0;
n10=0;
n13=0;
n14=0;
n15=0;
n16=0;
NZ=0;
DZ=0;
BZ=0;
OR=0;
NZero=0;
NNega=0;
None=0;
probability=zeros(NT(1,1),m);

%%
% "NZ": number of zero elements in a vector "s" and their indices.

% "DZ": number of non-zero elements in a vector "s" and their indices.

% "BZ": number of non-zero elements in a vector "s", and non-zero elements
% are not in the group related to action "c_i=a" and their indices.

% "OR": number of zero elements in a vector "s", and zero elements are not
% in the group related to action "c_i=a" and their indices.


while n6 <= m
    if s(n6)==0
        NZ=NZ+1;
        n7=n7+1;
        inNZ(1,n7)=n6;
    end
    if s(n6)~=0 
        n8=n8+1;
        DZ=DZ+1;
        inDZ(1,n8)=n6;
    end
    if s(n6)~=0 & n6~= a
        n9=n9+1;
        BZ=BZ+1;
        inBZ(1,n9)=n6;
    end
    if s(n6)==0 & n6~= a
        n10=n10+1;
        OR=OR+1;
        inOR(1,n10)=n6;
    end
    if s(a)==0 & m==2 & n6~= a
        n10=n10+1;
        OR=OR+1;
        inOR(1,n10)=n6;
    end
    n6=n6+1;
end

while n1<= NT(1,1)
    val=transition(n1,1:end)-s;
    while n2 <= m 
        if val(n2)==0 & type~=6
            NZero=NZero+1;
            n3=n3+1;
            inNZero(1,n3)=n2;
        end
        
        if type==6 
            if val(n2)==0 & n2~=a
                NZero=NZero+1;
                n3=n3+1;
                inNZero(1,n3)=n2;
            end
            if val(n2)==0 & n2==a
                NZero=NZero+1;
                n3=n3+1;
                inNZero(1,n3)=n2;
            end
        end
        if val(n2)==-1
            n4=n4+1;
            NNega=NNega+1;
            inNNega(1,n4)=n2;
        end
        if val(n2)==1
            n5=n5+1;
            None=None+1;
            inNone(1,n5)=n2;
        end
        n2=n2+1;
    end
    
    if type==1
        nothing(a,1)=1-lambda(a,1)-mu(a,1);
        nothing(inBZ,1)=1-mu(inBZ,1);
    end
    
    if type==2
        nothing(a,1)=1-lambda(a,1)-mu(a,1);
        nothing(inBZ,1)=1-mu(inBZ,1);
        nothing(inNZ,1)=1;
    end
    
    if type==3
        if BZ~=0
            nothing(a,1)=1-lambda(a,1)-mu(a,1);
            nothing(inBZ,1)=1-mu(inBZ,1);
            nothing(inNZ,1)=1;
        else
            nothing(a,1)=1-lambda(a,1)-mu(a,1);
            nothing(inNZ,1)=1;
        end
    end
    
    if type==4
        nothing(a,1)=1-lambda(a,1);
        nothing(inOR,1)=1;
    end
    if type==5
        if OR~=0 & m>= 3
            nothing(a,1)=1-lambda(a,1);
            nothing(inBZ,1)=1-mu(inBZ,1);
            nothing(inOR,1)=1;
        else
            nothing(a,1)=1-lambda(a,1);
            nothing(inBZ,1)=1-mu(inBZ,1);
        end
    end
    
    if type==6
        nothing(a,1)=1-mu(a,1)-lambda(a,1);
        nothing(inBZ,1)=1-mu(inBZ,1);
        nothing1(a,1)=1-mu(a,1);
        nothing1(inBZ,1)=1-mu(inBZ,1);
    end
    
    if type==7
        nothing(a,1)=1-mu(a,1)-lambda(a,1);
        nothing(inBZ,1)=1-mu(inBZ,1);
        nothing(inNZ,1)=1;
        nothing1(a,1)=1-mu(a,1);
        nothing1(inBZ,1)=1-mu(inBZ,1);
        nothing1(inNZ,1)=1;
    end
    
    if type==8
        nothing(a,1)=1-mu(a,1)-lambda(a,1);
        nothing(inNZ,1)=1;
        nothing1(a,1)=1-mu(a,1);
        nothing1(inNZ,1)=1;
    end
    
    if type==9
        nothing(a,1)=1-lambda(a,1);
        nothing(inBZ,1)=1-mu(inBZ,1);
        nothing(inNZ,1)=1;
        nothing1(a,1)=1-lambda(a,1);
        nothing1(inBZ,1)=1-mu(inBZ,1);
        if OR~=0
            nothing1(inOR,1)=1;
        end
        nothing2(a,1)=1;
        nothing2(inBZ,1)=1-mu(inBZ,1);
        nothing2(inNZ,1)=1;
    end
    
    if type==10
        nothing1(a,1)=1-lambda(a,1);
        nothing1(inBZ,1)=1-mu(inBZ,1);
        if OR~=0
            nothing1(inOR,1)=1;
        end
        nothing2(a,1)=1;
        nothing2(inBZ,1)=1-mu(inBZ,1);
        nothing2(inNZ,1)=1;
    end
    
    if type==11
        nothing(a,1)=1-mu(a,1);
        if OR~=0
            nothing(inOR,1)=1;
        end
        nothing2(a,1)=1;
        nothing2(inNZ,1)=1;
    end
    
     if None~=0 
         probability(n1,[inNone])=lambda(a,1);
     end
     
     if NZero~=0 & type<=5
         probability(n1,[inNZero])=nothing([inNZero],1);
     end
     
     if type==6 & NZero~=0 
         if val==0
             probability(n1,a)=nothing1(a,1);
             probability(n1,[inNZero])=nothing1([inNZero],1);
             
         elseif NZero~=0 & val(a)==0 
             probability(n1,a)=nothing(a,1); 
             probability(n1,[inNZero])=nothing([inNZero],1);
         elseif NZero~=0 & val(a)~=0
             probability(n1,[inNZero])=nothing([inNZero],1);
         end

     end
     
     if (type==7 | type==8) & NZero~=0
         if val==0
             probability(n1,a)=nothing1(a,1);
             probability(n1,[inNZero])=nothing1([inNZero],1);
         elseif NZero~=0 & val(a)==0
             probability(n1,a)=nothing(a,1); 
             probability(n1,[inNZero])=nothing([inNZero],1);
         elseif NZero~=0 & val(a)~=0
             probability(n1,[inNZero])=nothing([inNZero],1);
         end
     end
     
     if (type==9 | type==10) & NZero~=0
         if val==0
             probability(n1,a)=nothing2(a,1);
             probability(n1,[inNZero])=nothing2([inNZero],1);
         elseif NZero~=0 & val(a)==0 
             probability(n1,a)=nothing1(a,1); 
             probability(n1,[inNZero])=nothing1([inNZero],1);
         elseif NZero~=0 & val(a)~=0
             probability(n1,[inNZero])=nothing1([inNZero],1);
         end
     end
     
     if type==11
         if val==0
             probability(n1,a)=nothing(a,1);
             probability(n1,[inNZero])=nothing([inNZero],1);
         else
             probability(n1,[inNZero])=nothing([inNZero],1);
         end
     end
     if NNega~=0
         probability(n1,[inNNega])=mu([inNNega],1);
     end
    n1=n1+1;
    n2=1;
    n3=0;
    n4=0;
    n5=0;
    None=0;
    NZero=0;
    NNega=0;
    
    clear inNone inNNega inNZero 
end
NNega;
n13=0;
probability;
summation=0;
% checking that the summation of probabilities are 1!
while n13 < NT(1,1)
    n13=n13+1;
    summation=prod(probability(n13,1:end))+summation;
end
summation;
combined_probability=probability;

