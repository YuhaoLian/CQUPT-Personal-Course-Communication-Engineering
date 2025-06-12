clc;
clear;

%% ---------------- Queue Control Parameters ---------------%%
% Sim_Customer_Num=input('请输入仿真顾客总数'); %仿真顾客总数；
% Lambda=input('请输入仿真参数Lambda'); %到达模式服从参数为λ的泊松分布；
% Mu=input('请输入仿真参数Mu'); %服务时间服从参数为μ的指数分布；
%λ小于μ
Sim_Customer_Num = 10000;%仿真顾客总数；
Lambda=0.2; %到达模式服从参数为λ的泊松分布；
Mu=0.8;%服务时间服从参数为μ的指数分布；

%% ---------------- Initialization of Parameters ---------------%%
Time_Arrive=zeros(1,Sim_Customer_Num);
Time_Leave=zeros(1,Sim_Customer_Num);
Arrive_Num=zeros(1,Sim_Customer_Num);
Leave_Num=zeros(1,Sim_Customer_Num);
Interval_Arrive=zeros(1,Sim_Customer_Num);
Interval_Serve=zeros(1,Sim_Customer_Num);

Quene_Para = struct('Sim_Customer_Num', Sim_Customer_Num, 'Lambda', Lambda,'Mu',Mu, 'Time_Arrive', Time_Arrive,...
    'Time_Leave', Time_Leave, 'Arrive_Num', Arrive_Num, 'Leave_Num', Leave_Num, 'Interval_Arrive', Interval_Arrive,...
    'Interval_Serve', Interval_Serve);

save('Quene_Para', 'Quene_Para');


%% ---------------- Queue Generation ---------------%%
[t_Wait_avg,t_Queue_avg,QueLength_avg] = MM1_Queue(Quene_Para);

%% ---------------- μ Simulation ---------------%%
%---------------- Initialization of Simulation Parameters ---------------%
iter = 100;
Quene_Para.Lambda=0.2;
Mu_Simulation = 0.3:0.01:0.8;
Mu_Wait_avg_iter = zeros(iter,size(Mu_Simulation,2));
Mu_Queue_avg_iter = zeros(iter,size(Mu_Simulation,2));
Mu_QueLength_avg_iter = zeros(iter,size(Mu_Simulation,2));

for i = 1:iter
    Mu_Wait_avg = zeros(size(Mu_Simulation));
    Mu_Queue_avg = zeros(size(Mu_Simulation));
    Mu_QueLength_avg = zeros(size(Mu_Simulation));
    for j =1:size(Mu_Simulation,2)
        fprintf('== Mu is %g and iteration is %g == \n',Mu_Simulation(j),i)
        Quene_Para.Mu=Mu_Simulation(j);
        [Mu_Wait_avg(j),Mu_Queue_avg(j),Mu_QueLength_avg(j)] = MM1_Queue(Quene_Para);
    end
    Mu_Wait_avg_iter(i,:)=Mu_Wait_avg;
    Mu_Queue_avg_iter(i,:)=Mu_Queue_avg;
    Mu_QueLength_avg_iter(i,:)=Mu_QueLength_avg;
end
Mu_Wait_avg = sum(Mu_Wait_avg_iter)/iter;
Mu_Queue_avg = sum(Mu_Queue_avg_iter)/iter;
Mu_QueLength_avg = sum(Mu_QueLength_avg_iter)/iter;

figure;
subplot(1,3,1);
plot(Mu_Simulation,Mu_Wait_avg,'LineWidth',1.5,'MarkerSize',12);
xlabel('\mu');
ylabel('平均逗留时间(t/h)');
title('仿真平均逗留时间');
subplot(1,3,2);
plot(Mu_Simulation,Mu_Queue_avg,'LineWidth',1.5,'MarkerSize',12);
xlabel('\mu');
ylabel('平均排队时间(t/h)');
title('仿真平均排队时间');
subplot(1,3,3);
plot(Mu_Simulation,Mu_QueLength_avg,'LineWidth',1.5,'MarkerSize',12);
xlabel('\mu');
ylabel('平均队长(t/h)');
title('仿真系统中平均队长');

%% ---------------- λ Simulation ---------------%%
%---------------- Initialization of Simulation Parameters ---------------%
iter = 100;
Quene_Para.Mu=0.8;
Lambda_Simulation = 0.1:0.01:0.6;
Lambda_Wait_avg_iter = zeros(iter,size(Lambda_Simulation,2));
Lambda_Queue_avg_iter = zeros(iter,size(Lambda_Simulation,2));
Lambda_QueLength_avg_iter = zeros(iter,size(Lambda_Simulation,2));

for i = 1:iter
    Lambda_Wait_avg = zeros(size(Lambda_Simulation));
    Lambda_Queue_avg = zeros(size(Lambda_Simulation));
    Lambda_QueLength_avg = zeros(size(Lambda_Simulation));
    for j =1:size(Lambda_Simulation,2)
        fprintf('== Lambda is %g and iteration is %g == \n',Lambda_Simulation(j),i)
        Quene_Para.Lambda=Lambda_Simulation(j);
        [Lambda_Wait_avg(j),Lambda_Queue_avg(j),Lambda_QueLength_avg(j)] = MM1_Queue(Quene_Para);
    end
    Lambda_Wait_avg_iter(i,:)=Lambda_Wait_avg;
    Lambda_Queue_avg_iter(i,:)=Lambda_Queue_avg;
    Lambda_QueLength_avg_iter(i,:)=Lambda_QueLength_avg;
end
Lambda_Wait_avg = sum(Lambda_Wait_avg_iter)/iter;
Lambda_Queue_avg = sum(Lambda_Queue_avg_iter)/iter;
Lambda_QueLength_avg = sum(Lambda_QueLength_avg_iter)/iter;

figure;
subplot(1,3,1);
plot(Lambda_Simulation,Lambda_Wait_avg,'r','LineWidth',1.5,'MarkerSize',12);
xlabel('\lambda');
ylabel('平均逗留时间(t/h)');
title('仿真平均逗留时间');
subplot(1,3,2);
plot(Lambda_Simulation,Lambda_Queue_avg,'r','LineWidth',1.5,'MarkerSize',12);
xlabel('\lambda');
ylabel('平均排队时间(t/h)');
title('仿真平均排队时间');
subplot(1,3,3);
plot(Lambda_Simulation,Lambda_QueLength_avg,'r','LineWidth',1.5,'MarkerSize',12);
xlabel('\lambda');
ylabel('平均队长(t/h)');
title('仿真系统中平均队长');