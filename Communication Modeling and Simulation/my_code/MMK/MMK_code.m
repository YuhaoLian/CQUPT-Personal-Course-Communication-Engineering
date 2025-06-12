clc;
clear;

%% ---------------- Queue Control Parameters ---------------%%
% Sim_Customer_Num=input('请输入仿真顾客总数'); %仿真顾客总数；
% Lambda=input('请输入仿真参数Lambda'); %到达模式服从参数为λ的泊松分布；
% Mu=input('请输入仿真参数Mu'); %服务时间服从参数为μ的指数分布；
%λ小于μ
Sim_Customer_Num = 10000000;%仿真顾客总数；
Lambda=0.9; %到达模式服从参数为λ的泊松分布；
Mu=0.4;%服务时间服从参数为μ的指数分布；
K=3;%服务台数量K

%% ---------------- Initialization of Parameters ---------------%%
Time_Arrive=zeros(1,Sim_Customer_Num);
Time_Leave=zeros(1,Sim_Customer_Num);
Interval_Arrive=zeros(1,Sim_Customer_Num);
Interval_Serve=zeros(1,Sim_Customer_Num);

Quene_Para = struct('Sim_Customer_Num', Sim_Customer_Num, 'Lambda', Lambda,'Mu',Mu,'K',K, 'Time_Arrive', Time_Arrive,...
    'Time_Leave', Time_Leave, 'Interval_Arrive', Interval_Arrive,'Interval_Serve', Interval_Serve);

save('Quene_Para', 'Quene_Para');


%% ---------------- Queue Generation ---------------%%
%---------------- MM3 Queue ---------------%
[t_Wait_avg,t_Queue_avg,QueLength_avg] = MMK_Queue(Quene_Para);

%---------------- 3*MM1 Queue ---------------%
Quene_Para.K = 1;
Quene_Para.Lambda = 0.3;
[t_Wait_avg_1,t_Queue_avg_1,QueLength_avg_1] = MMK_Queue(Quene_Para);
