
%%================ Exact DP

clear all
clc
tic

m=input('Enter the number of prices:')
N=input('Enter the number of resources:')



%% Time Horizon, price vector,
% T: variable for the horizon
%  C: Price vector

load lambda.mat
load mu.mat
%load T.mat
load C.mat

T=200;

%lambda=[0.4;0.5;0.4;0.3]
%mu=[0.35;0.4;0.45;0.25]
%C=[0.9 1 1.1 1.15]
%e.g. lambda=[0.6;0.5;0.3]; "m*1" vector
%e.g. mu=[0.2;0.2;0.4]; "m*1" vector

%e.g. T=60;
%e.g. C=[0.9 1 1.1]; "1*m" vector

%%
lambda = [0.090596; 0.048632; 0.015657; 0.005088]
mu = [0.483723; 0.444019; 0.024843; 0.335103]
C = [0.8 0.9 1 1.1]

%Suitable lambda and mu found:
%lambda=[0.049297; 0.071014; 0.084126; 0.098124]
%mu= [0.158740; 0.158214; 0.108782; 0.125521]
%State with maximum difference: [0, 0, 0, 0]
%Maximum difference in transition probabilities: 0.0981


%%
%  Maximum difference in transition probabilities: 0.0982

%lambda = [0.026091; 0.007525; 0.098181]
%mu = [0.263856; 0.073836; 0.036567]
%C = [0.8 0.9 1]

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

%% In this subsection we find the value of the states at the terminal period
% for all the states

% Addition(s,C): It is a written function whose inputs are a state "s" and
% the price vector "C" e.g. Addition([1 2],[0.9 1.1]), which calculates the
% value of a state.

% Terminal_revenue: It is a variable to save the value function of each
% state at the terminal time step, its size is "NS(1,1)*(m+1)" ;
% each row has the e.g. [s,s*C'], "C'" denotes transpose of vector price
% "C".

n1=1;
while n1<= NS(1,1)
    Addition(S(n1,1:end),C);
    Terminal_revenue(n1,1:(m+1),1)=[S(n1,1:end) Addition(S(n1,1:end),C)];
    n1=n1+1;
    
end


%% In this subsection we apply DP for each state and the entire horizon

% we use variable "a" as an index variable indicating a price "c_i".
% "k" is the time slot counter, k={0,1,...,T}.
% "Revenue" is the variable that saves the value function for each time
% slot "k", it is a 3 dimension matrix; also it saves the optimal decision 
% for each state 
% e.g. Revenue(1:NS,1:m+2,k), the row index points at a state, each row
% contains the corresponding state, Decision and optimal revenue,
% respectively.
% 

Discount_factor=0.9;

a=1;
n1=1;
n2=1;
k=T-1;
Revenue(1:NS(1,1),1:(m+1),T)=Terminal_revenue(1:NS(1,1),1:(m+1));
Revenue(1:NS(1,1),(m+2),T)=0;

%=============================

% This subsection is the main loop of DP
% we check whether or not each price is zero 
% loop on the variable "k={0,1,...,T}" represents the horizon which is the backward
% procedure; loop on variable "n1" denotes the index of state space that evaluate
% the policy for each state; loop on variable "a={1,2,...,m}" dentoes evaluation "m" 
% different prices for each state

%=====================

% stateanalysis(s,N,m,a): It is a function whose inputs are the state "s",
% the number of resources "N", number of prices "m", and an action "a=c_i"; whose
% output are possible transitions from a state if an action "a=c_i" has been taken;
% the other output is the probability of going from a state to other state by 
% taking an action "a=c_i", 
% e.g. [transition combined_probability]=stateanalysis([1 1],3,2,1)

%=====================
% 
% In DP we save the result as a look up table in the matrix "Revenue"
% "pair_finding": given a state at time slot "k", having possible
% "transitions", and "combined_probability", this function helps to find the 
% "Revenue" at the time slot "k+1" corersponding to each state "x(k+1)"
% based on "transition" matrix from matrix "Revenue".
% 
% possible_Revenue(1,a): It is a vector that saves the value fucntion for
% making an action "c_i=a" at each state, its size is 1*m i.e. at each state we evaluate "m"
% different policies and save them into this variable, the DP chooses the best
% policy from this vector for each state at each time step and inserts
% into the variable "Revenue"

%================================================

while k>= 1
    
    %lambda=[0.5+0.1*cos(k*pi/3);0.3+.1*sin(k*pi/7);0.2+0.1*cos(k*pi/5)];
    %mu=[0.3-0.1*cos(k*pi/3);0.4-0.1*sin(k*pi/7);0.5-0.1*cos(k*pi/5)];

    %lambda=[0.5+0.1*cos(k*pi/3);0.3+.1*sin(k*pi/7);0.2+0.1*cos(k*pi/5); 0.15-0.05*cos(k*pi/6)];
    %mu=[0.3-0.1*cos(k*pi/3);0.4-0.1*sin(k*pi/7);0.5-0.1*cos(k*pi/5); 0.6+0.05*cos(k*pi/6)];
    
    while n1<= NS(1,1) 
        while a <= m 
            if C(a)~=0
                S(n1,1:end);
                [transition combined_probability]=stateanalysis(S(n1,1:end),N,m,a,lambda,mu);
                x=Revenue(1:end,1:end,k+1);
                %Revenue(1:end,1:end,k+1)
                Rev=Pair_finding(transition,x,NS(1,1),combined_probability);
                possible_Revenue(1,a)=Rev;
            end
            a=a+1;
        end
            
            [max_revenue_state price]=max(possible_Revenue);
            %price = randi([1 4],1,1);

            %max_revenue_state = possible_Revenue(price);

            Revenue(n1,(m+1),k)=Terminal_revenue(n1,(m+1))+Discount_factor*max_revenue_state;
            Revenue(n1,1:m,k)=S(n1,1:end);
            Revenue(n1,(m+2),k)=price;
        n1=n1+1;
        a=1;
    end
    k=k-1;
    n1=1;
end

%%

%% Frequency Distribution Before Reservation

%m=4;
st_size(1)=NS(1);
T=200;
for n40=1:T
    accumulative_var=zeros(1,m);
    for a=1:m
        for i=1:st_size(1)
            if Revenue(i,m+2,n40)==a
                accumulative_var(1,a)=accumulative_var(1,a)+1;
                frequency_decision_BD(n40,a)=accumulative_var(1,a);
            end 
        end
    end
    %accumulative_var=zeros(1,m+2);
end

%%
hold on

subplot(3,1,1)
%plot(frequency_decision_net(:,1)/st_size(1))
hold on
plot(frequency_decision_BD(:,1)/st_size(1))
legend('Booking Case','No Booking Case')
subplot(3,1,2)
%plot(frequency_decision_net(:,2)/st_size(1))
hold on
plot(frequency_decision_BD(:,2)/st_size(1))
subplot(3,1,3)
%plot(frequency_decision_net(:,3)/st_size(1))
hold on
plot(frequency_decision_BD(:,3)/st_size(1))
subplot(4,1,4)
%plot(frequency_decision_net(:,4)/st_size(1))
%hold on
%plot(frequency_decision_BD(:,4)/st_size(1))
%%
% Set the figure size to be larger
figure('Position', [100, 100, 1000, 800]);

% First subplot
subplot(3,1,1)
hold on
plot(frequency_decision_BD(:,1)/st_size(1), 'LineWidth', 3)
%plot(frequency_decision_net(:,1)/st_size(1), 'LineWidth', 3) % Uncomment if you want to include this
ylabel('C_1', 'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold')
%legend('Booking Case', 'No Booking Case', 'Location', 'northeast', 'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold')
grid on

% Second subplot
subplot(3,1,2)
hold on
plot(frequency_decision_BD(:,2)/st_size(1), 'LineWidth', 3)
%plot(frequency_decision_net(:,2)/st_size(1), 'LineWidth', 3) % Uncomment if you want to include this
ylabel('C_2', 'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold')
grid on

% Third subplot
subplot(3,1,3)
hold on
plot(frequency_decision_BD(:,3)/st_size(1), 'LineWidth', 3)
%plot(frequency_decision_net(:,3)/st_size(1), 'LineWidth', 3) % Uncomment if you want to include this
xlabel('Time Slot', 'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold')
ylabel('C_3', 'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold')
grid on

% Optional Fourth subplot
% subplot(4,1,4)
% hold on
% plot(frequency_decision_BD(:,4)/st_size(1), 'LineWidth', 3)
% plot(frequency_decision_net(:,4)/st_size(1), 'LineWidth', 3) % Uncomment if you want to include this
% xlabel('Time Slot', 'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold')
% ylabel('Frequency', 'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold')
% grid on

% General settings
set(gcf, 'Color', 'w'); % Set the background to white

% Apply font settings to all axes
set(findall(gcf, 'Type', 'axes'), 'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold');


%%

% Assuming Revenue is already defined and has size [56, 5, 200]

% Extract the data you want to plot
data_to_plot = squeeze(Revenue(12, 5, :));

% Plotting the data
figure;
plot(data_to_plot);
xlabel('Time Step');
ylabel('Value');
title('Plot of Revenue(1, 5, :)');



%%

toc

