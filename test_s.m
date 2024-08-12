
clear all


% Given values
lambda = [0.3; 0.4; 0.2]; % Lambda values
mu = [0.2; 0.3; 0.1]; % Mu values

% Given values
lambda = [0.2; 0.2; 0.2]; % Lambda values
mu = [0.2; 0.2; 0.2]; % Mu values


m = 3
N= 3

%%

%% State Space Construction
% S: state_space
% Here we construct a matrix S whose rows contain all the possible combinations 
% of at most N objects on m places.
% To that aim, we use a routine base(m,N,n) which computes the representation of number
% n in the base N+1 over m figures, then we cancel the combinations of
% figures whose sum exceeds N
% 
n=0;
l=0;
while n<=((N+1)^m - 1)
    s=base(m,N,n);
    if sum(s)<= N
    l=l+1;
    S(l,1:m)=s;
    end
    n=n+1;
end

% NS: size of entire state space
NS=size(S);


%%



% Compute the maximum difference in transition probabilities
max_diff = max_probability_difference_all(S, N, m, lambda, mu);


