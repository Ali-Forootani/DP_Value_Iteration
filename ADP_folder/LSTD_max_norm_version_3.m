%%
clear all
clc

l=0;
m=input('Enter the number of prices:')
N=input('Enter the parking capacity:')

%initial_state=[1 2 1 1];
%price_vector=[0.9 1 1.1];


initial_state=[0 0 0];

s=initial_state;
a=1;
n=0;
%path(n,1:m)=s;
dis_factor=0.9;
horizon=1;
summation_c=0;
summation_d=0;
Decision=[];

%lambda=[0.55;0.5;0.3];
%mu=[0.15;0.2;0.35];

lambda=[0.6;0.5;0.3];
mu=[0.2;0.2;0.4];

price_vector=[0.9 1 1.1];


%%
lambda = [0.090596; 0.048632; 0.015657; 0.005088]
mu = [0.483723; 0.444019; 0.024843; 0.335103]

price_vector = [0.8 0.9 1 1.1]


%r_hat =[725; 3.84; 5.02; 2.55];

r_hat=[1; 0.8; 0.9; 1; 1.1];

r_max_col=[];
%r_bar=montecarlo_r_hat(m,N,initial_state,price_vector,horizon,dis_factor);

k=1;
%%
l_m=2;

while l_m<= 50000

sd=0;

mod5=0;


while sd==0
    
    if l_m == 2

    rand_num=randi((N+1)^m-1,1,1);
    possible_initial_state=base(m,N,rand_num);
    if sum(possible_initial_state)<= N
        sd=1;
        initial_state=possible_initial_state;
    end

    end

    if l_m > 2
        a_l = randi(m, 1, 1);
        [trans_initial_state combined_probability] = stateanalysis(s, N, m, a_l, lambda, mu);
        s_trans_init = size(trans_initial_state);
        a_d = randi(s_trans_init(1), 1, 1);
        initial_state = trans_initial_state(a_d,:);
        sd=1;

    end

end



n20=1;
s=initial_state;

indicator_inequality=0;

%path(n20,1:m)=s;

if l_m==1
    r_bar=montecarlo_r_hat(m,N,initial_state,price_vector,horizon,dis_factor);
end

while k<= horizon
    %%    
    f=rand;
    s;
    
   

    
    %sum(s)
    main_path(k,1:m)=s;
    
    if sum(s)< N
        
        while a<=m
    %%      
            s;
            n=n+1;
            a;
            f;
            [transition combined_probability]=stateanalysis(s,N,m,a,lambda,mu);
           
            [next_state]=path_creation(transition,combined_probability,f,a);
            possible_nextstate(a,:)=next_state;
            
            
            %%
            %Function to calculate "Phi" according to linear approximation of value
            % function
            ph_k=phi_cal2(m,s);
            ph_kprime=phi_cal2(m,next_state);
    
            %%
            [c_k summation_ck]=C_cal(ph_k,ph_kprime,dis_factor,k,summation_c);
            summation_c=summation_ck;
            %%
            g_state=g_cal(s,price_vector);
            [d_k summation_dk]=d_cal(ph_k,k,summation_d,s,g_state);
            summation_d=summation_dk;
            %%
            Sigma=2*eye(length(c_k));
            beta=0.01;
            %beta=2;
            bias=beta*eye(length(c_k));
            if l_m==1
                r_hat=r_bar*rand(length(c_k),1);
            end
            r(:,a)=inv(c_k'*inv(Sigma)*c_k+bias)*(c_k'*inv(Sigma)*d_k+bias*r_hat);
            r_norm(a)=norm(r(:,a));
            
            
            %%   matrix c 
            
            cc_k(:,:,a)=c_k;
            
            
            %%
            s;
            
            a=a+1;
            indicator_rej=1;
        end
        
%         if indicator_rej==1
%             n=n+1;
%             [transition_re combined_probability_re]=stateanalysis_re(s,N,m,mu);
%             %g=size(s);
%             s;
%             [next_state_re]=path_creation_re(transition_re,combined_probability_re,f);
%             possible_nextstate(m+1,:)=next_state_re;
%             
%             ph_k=phi_cal2(m,s);
%             ph_kprime=phi_cal2(m,next_state_re);
%             
%             [c_k summation_ck]=C_cal(ph_k,ph_kprime,dis_factor,k,summation_c);
%             summation_c=summation_ck;
%             
%             g_state=g_cal(s,price_vector);
%             [d_k summation_dk]=d_cal(ph_k,k,summation_d,s,g_state);
%             summation_d=summation_dk;
%             
%             Sigma=2*eye(length(c_k));
%             beta=0.01;
%             bias=beta*eye(length(c_k));
%             if l_m==1
%                 r_hat=r_bar*rand(length(c_k),1);
%             end
%             
%             r(:,m+1)=inv(c_k'*inv(Sigma)*c_k+bias)*(c_k'*inv(Sigma)*d_k+bias*r_hat);
%             r_norm(m+1)=norm(r(:,m+1));
%             
%             path(n,1:m)=s;
%             
%         end
%         
        k_normal=k;
        path(k,1:m,l_m)=s;
        
        for i=1:m
                instant_rev(i)=g_state+dis_factor*phi_cal2(m,possible_nextstate(i,:))'*r(:,i);
        end
        
        
        current_state=s;
        r_norm_path(k,1:m,l_m)=r_norm;
        
        %% Control selection
        %
        %Max
        [r_max_instant(k) Decision(k)]=max(instant_rev);
        Decision(k);
        
        %Min
        %[r_max_instant(k) Decision(k)]=min(instant_rev);
        
        %Random
        %random_deci=randi([1 m],1,1);
        %r_norm_max(k)=instant_rev(:,random_deci);
        %Decision(k)=random_deci;
        
        %Something else
          % f_x=rand;
          % if f_x < 0.6
          %     [min_death associate]=min(mu);
          %     Decision(k)=associate;
          %     r_max_instant(k)=instant_rev(:,associate);
          % else
          %     random_deci=randi([1 m],1,1);
          %     r_norm_max(k)=instant_rev(:,random_deci);
          %     Decision(k)=random_deci;
          % end
%        
         %The other policy
%         Decision(k)=mod(sum(s),3)+1;
%         r_norm_max(k)=instant_rev(:,Decision(k));
         
        %%
        s;
        r_max(:,k)=r(:,Decision(k));
        next_state=possible_nextstate(Decision(k),:);
        
        %About matrix c_k
        
        c_k_max(:,:,k)=cc_k(:,:,Decision(k));
        c_k_max_eig(:,k)=eig(c_k_max(:,:,k));
        %
        
        s=next_state;
        main_path_iteration(k,1:m+1,l_m)=[current_state Decision(k)];
        
        %==========
        %main_path(k,1:m)=path(Decision(k),1:m,l_m);
        %========
        
        %main_path_iteration(k,:,l_m)=path(Decision(k),1:m);
        r_hat=r_max(:,k);
        
        indicator_inequality=1;
        
        main_decision(k,l_m)=Decision(k);
    end
    
%%
if sum(s)==N & indicator_inequality==0
            s;
            %n=n+1;
            %a;
            f;
            [transition_ca combined_probability_ca]=stateanalysis_re(s,N,m,mu);
            %g=size(s);
 
            [next_state_ca]=path_creation_re(transition_ca,combined_probability_ca,f);
            %possible_nextstate(a,:)=next_state_ca;
            %%
            ph_k=phi_cal2(m,s);
            ph_kprime=phi_cal2(m,next_state_ca);
            
            %%
            [c_k summation_ck]=C_cal(ph_k,ph_kprime,dis_factor,k,summation_c);
            summation_c=summation_ck;
            %%
            g_state=g_cal(s,price_vector);
            [d_k summation_dk]=d_cal(ph_k,k,summation_d,s,g_state);
            summation_d=summation_dk;
            
            %%
            Sigma=2*eye(length(c_k));
            beta=0.01;
            bias=beta*eye(length(c_k));
            if l_m==1
                r_hat=r_bar*rand(length(c_k),1);
            end
            r_ca(:,1)=inv(c_k'*inv(Sigma)*c_k+bias)*(c_k'*inv(Sigma)*d_k+bias*r_hat);
            r_norm_ca(1)=norm(r_ca(:,1));
            
            k_un=k;
            path_ca(k,1:m,l_m)=s;
            
            
            
            %%
            current_state_ca=s;
            r_norm_path(k,1:m,l_m)=r_norm_ca;
            %next_state=next_state_ca;
            Decision_ca(k)=m+1;
            %r_max(:,k)=r_ca(:,1);
            s=next_state_ca;
            main_path_iteration(k,:,l_m)=[current_state_ca Decision_ca(k)];
            
            %============
            %main_path(k,1:m)=path_ca(1,1:m);
            %============
            
            r_hat=r_ca(:,1);
            r_max(:,k)=r_ca(:,1);
            
            main_decision(k,l_m)=Decision_ca(k);
            
            %About matrix c_k
            c_k_max(:,:,k)=c_k;
            c_k_max_eig(:,k)=eig(c_k);
end


% 
%status=[current_state;next_state];

%main_path(k,1:m)=current_state;


%current_state;

%Decision(k);


indicator_inequality=0;


a=1;
n=0;
k=k+1;
indicator_rej=0;

%path=[];
%path_ca=[];

possible_nextstate=[];

% current_state_ca=[];
% next_state_ca=[];
% next_state=[];
% current_state=[];

%fprintf('*********');

end

%main_path=[];

r_max_col=[r_max_col r_max];

Decision=[];
Decision_ca=[];

summation_c;


c_k_max_eig_traj(:,:,l_m)=c_k_max_eig(:,:);

c_k_max_eig=[];


summation_d;

k=1;

l_m=l_m+1;

%path=[];

%path_ca=[];

possible_nextstate=[];

end
%stateanalysis()


r_s_m=size(r_max_col);
n50=0;
n70=1;
while n70<=r_s_m(2)
    if norm(r_max_col(:,n70)) > 1
        n50=n50+1;
        r_max_f(n50)=norm(r_max_col(:,n70));
    end
    n70=n70+1;
end

%%

c_k_max_eig_traj=real(c_k_max_eig_traj);

interme_traj=size(c_k_max_eig_traj);

length_traj=interme_traj(2)*interme_traj(3);
n75=1;
n76=1;
n77=1;
while n75<= interme_traj(2)
    
    while n76<= interme_traj(3)
    matrix_c_traj(:,n77) = c_k_max_eig_traj(:,n75,n76);
    n76=n76+1;
    n77=n77+1;
    end
    n75=n75+1;
    n76=1;
end


