function [t_Wait_avg,t_Queue_avg,QueLength_avg] = MM1_Queue(Quene_Para)

Sim_Customer_Num = Quene_Para.Sim_Customer_Num;
Lambda=Quene_Para.Lambda; 
Mu=Quene_Para.Mu;
Time_Arrive=Quene_Para.Time_Arrive;
Time_Leave=Quene_Para.Time_Leave;
Arrive_Num=Quene_Para.Arrive_Num;
Leave_Num=Quene_Para.Leave_Num;
Interval_Arrive=Quene_Para.Interval_Arrive;
Interval_Serve=Quene_Para.Interval_Serve;

Interval_Arrive=-log(rand(1,Sim_Customer_Num))/Lambda;%到达时间间隔
Interval_Serve=-log(rand(1,Sim_Customer_Num))/Mu;%服务时间
Time_Arrive(1)=Interval_Arrive(1);%顾客到达时间

ArriveNum(1)=1;
for i=2:Sim_Customer_Num
    Time_Arrive(i)=Time_Arrive(i-1)+Interval_Arrive(i);
    ArriveNum(i)=i;
end

Time_Leave(1)=Time_Arrive(1)+Interval_Serve(1);%顾客离开时间
LeaveNum(1)=1;
for i=2:Sim_Customer_Num
    if Time_Leave(i-1)<Time_Arrive(i)
        Time_Leave(i)=Time_Arrive(i)+Interval_Serve(i);
    else
        Time_Leave(i)=Time_Leave(i-1)+Interval_Serve(i);
    end
    LeaveNum(i)=i;
end

t_Wait=Time_Leave-Time_Arrive; %各顾客在系统中的逗留时间
t_Wait_avg=mean(t_Wait);
t_Queue=t_Wait-Interval_Serve;%各顾客在系统中的排队时间
t_Queue_avg=mean(t_Queue);

Timepoint=[Time_Arrive,Time_Leave];%系统中离开事件、到达事件发生的时间
Timepoint=sort(Timepoint);
ArriveFlag=zeros(size(Timepoint));%到达时间标志
CusNum=zeros(size(Timepoint));%到达或者离开事件发生时系统的顾客数量
temp=2;
CusNum(1)=1;
for i=2:length(Timepoint)
    if (temp<=length(Time_Arrive))&&(Timepoint(i)==Time_Arrive(temp))%这里隐含着任何事件不可能同时发生的假设
        CusNum(i)=CusNum(i-1)+1;
        temp=temp+1;
        ArriveFlag(i)=1;
    else
        CusNum(i)=CusNum(i-1)-1;
    end
end
%计算任何相邻事件的时间间隔
Time_interval=zeros(size(Timepoint));
Time_interval(1)=Time_Arrive(1);
for i=2:length(Timepoint)
    Time_interval(i)=Timepoint(i)-Timepoint(i-1);
end
%系统中平均顾客数计算
CusNum_fromStart=[0 CusNum];
CusNum_avg=sum(CusNum_fromStart.*[Time_interval 0] )/Timepoint(end);

%计算事件发生时系统队列的长度
QueLength=zeros(size(CusNum));
for i=1:length(CusNum)
    if CusNum(i)>=2
        QueLength(i)=CusNum(i)-1;%系统中只有1个服务台
    else
        QueLength(i)=0;
    end
end
QueLength_avg=sum([0 QueLength].*[Time_interval 0] )/Timepoint(end);%系统平均队长

%仿真值与理论值比较
% disp(['理论平均逗留时间=',num2str(1/(Mu-Lambda))]);
% disp(['理论平均排队时间=',num2str(Lambda/(Mu*(Mu-Lambda)))]);
% disp(['理论系统中平均队长=',num2str(Lambda*Lambda/(Mu*(Mu-Lambda)))]);
% 
% disp(['仿真平均逗留时间=',num2str(t_Wait_avg)])
% disp(['仿真平均排队时间=',num2str(t_Queue_avg)])
% disp(['仿真系统中平均队长=',num2str(QueLength_avg)]);
end

