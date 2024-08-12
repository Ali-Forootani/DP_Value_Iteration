function [Next_Revenue] = Pair_finding(transition,x,NS,combined_probability)

% This function takes each state and its probability, corresponding to each row of
% "transition" and "combined_probability" matrices and finds the state
% inside the matrix "Revenue=x" along with its value function.

% matrix "x" has the elements of matrix "Revenue" at the time slot "k+1" 

NT=size(transition);
n1=1;
n2=1;
n3=0;
while n1<= NT(1,1)
    while n2<= NS
        b=isequal(transition(n1,1:NT(1,2)),x(n2,1:NT(1,2)));
        if b==1
            n3=n3+1;
            betta=x(n2,end-1);
            sum_Next_Revenue(1,n3)=x(n2,end-1)*prod(combined_probability(n1,1:NT(1,2)));
        end
        n2=n2+1;
    end
    
    n1=n1+1;
    n2=1;
end
Next_Revenue=sum(sum_Next_Revenue);