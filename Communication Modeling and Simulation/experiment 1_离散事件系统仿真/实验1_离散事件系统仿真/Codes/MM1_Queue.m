clear;
clc;
%M/M/1排队系统仿真
SimTotal=100000; %仿真顾客总数
Lambda=0.3; %到达率Lambda
Mu=0.4; %服务率Mu；
t_Arrive=zeros(1,SimTotal); 
t_Leave=zeros(1,SimTotal);
ArriveNum=zeros(1,SimTotal);
LeaveNum=zeros(1,SimTotal);

Interval_Arrive=-log(rand(1,SimTotal))/Lambda;%到达时间间隔
Interval_Serve=-log(rand(1,SimTotal))/Mu;%服务时间
t_Arrive(1)=Interval_Arrive(1);%顾客到达时间
ArriveNum(1)=1;
for i=2:SimTotal
    t_Arrive(i)=t_Arrive(i-1)+Interval_Arrive(i);
    ArriveNum(i)=i;
end

t_Leave(1)=t_Arrive(1)+Interval_Serve(1);%顾客离开时间
LeaveNum(1)=1;
for i=2:SimTotal
    if t_Leave(i-1)<t_Arrive(i)
        t_Leave(i)=t_Arrive(i)+Interval_Serve(i);
    else
        t_Leave(i)=t_Leave(i-1)+Interval_Serve(i);
    end
    LeaveNum(i)=i;
end

t_Wait=t_Leave-t_Arrive; %各顾客在系统中的逗留时间
t_Wait_avg=mean(t_Wait);
t_Queue=t_Wait-Interval_Serve;%各顾客在系统中的排队时间
t_Queue_avg=mean(t_Queue);

Timepoint=[t_Arrive,t_Leave];%系统中离开事件、到达事件发生的时间
Timepoint=sort(Timepoint);
ArriveFlag=zeros(size(Timepoint));%到达时间标志
CusNum=zeros(size(Timepoint));%到达或者离开事件发生时系统的顾客数量
temp=2;
CusNum(1)=1;
for i=2:length(Timepoint)
    if (temp<=length(t_Arrive))&&(Timepoint(i)==t_Arrive(temp))%这里隐含着任何事件不可能同时发生的假设
        CusNum(i)=CusNum(i-1)+1;
        temp=temp+1;
        ArriveFlag(i)=1;
    else
        CusNum(i)=CusNum(i-1)-1;
    end
end
%计算任何相邻事件的时间间隔
Time_interval=zeros(size(Timepoint));
Time_interval(1)=t_Arrive(1);
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
disp(['理论平均逗留时间=',num2str(1/(Mu-Lambda))]);
disp(['理论平均排队时间=',num2str(Lambda/(Mu*(Mu-Lambda)))]);
disp(['理论系统中平均队长=',num2str(Lambda*Lambda/(Mu*(Mu-Lambda)))]);

disp(['仿真平均逗留时间=',num2str(t_Wait_avg)])
disp(['仿真平均排队时间=',num2str(t_Queue_avg)])
disp(['仿真系统中平均队长=',num2str(QueLength_avg)]);
