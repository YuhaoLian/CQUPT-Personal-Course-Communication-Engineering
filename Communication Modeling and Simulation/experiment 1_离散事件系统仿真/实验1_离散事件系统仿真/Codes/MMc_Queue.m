clear;clc;close all;
% M/M/K排队系统仿真
% 基于事件列表模型

SimTotal=100000; %仿真顾客总数
Lambda=0.9; %到达率Lambda
Mu=0.4; %服务率Mu；
c = 3; %服务台数量

%按照分布，初始化到达间隔和服务间隔
Interval_Arrive=-log(rand(1,SimTotal))/Lambda;
Interval_Serve=-log(rand(1,SimTotal))/Mu;
%通过累计求和，计算顾客的达到时刻；
t_Arrive = cumsum(Interval_Arrive);

%初始化离开时刻
t_Leave=zeros(1,SimTotal);
%定义时间列表的三个参量的索引值：
%事件时刻，事件类型（1代表到达，-1代表离开），客户编号
ETIME = 1; ETYPE = 2; ECUST_NO = 3;
%初始化事件列表，加入第一个客户到达的事件
evList(1,ETIME) = t_Arrive(1);
evList(1,ETYPE) = 1; 
evList(1,ECUST_NO) = 1;
%初始化队列
qList = [];
%下一个到达的客户编号；
nextN = 2;
%繁忙服务台数量
busyN = 0;
%按照时间加权求和的队列长度
nQ = 0;
%队列长度上一次变化的时刻
tQ = 0;
%主逻辑，以事件列表非空作为循环条件
while(~isempty(evList))
    %到达事件
    if (evList(1,ETYPE) == 1)
        %如果有服务台空闲
        if (busyN < c)
            %繁忙服务台加一
            busyN = busyN + 1;
            %设置该客户的离开事件，按逻辑应该是删除到达事件，再添加离开事件
            %此处优化为直接修改事件类型和时刻，
            evList(1,ETIME ) = evList(1,ETIME) + Interval_Serve(evList(1,ECUST_NO));
            evList(1,ETYPE) = -1;
        %服务台没有空闲
        else
            %队列长度将发生变化，累积前一段时刻的队列长度求和
            nQ = nQ + length(qList)*(evList(1,ETIME) - tQ);
            tQ = evList(1,ETIME);
            %将顾客排入队列
            qList(end+1) = evList(1, ECUST_NO);
            %删除该顾客达到事件
            evList(1,:) = [];
        end
        %如果还有客户将要到达，处理新顾客到达事件，加入事件列表
        if (nextN <= SimTotal)
            evList(end+1, ETIME) = t_Arrive(nextN);
            evList(end, ETYPE) = 1;
            evList(end, ECUST_NO) = nextN;
            nextN = nextN + 1;
        end
    %离开事件
    else
        %繁忙服务台减一
        busyN = busyN - 1;
        %记录该顾客的离开事件
        t_Leave(evList(1, ECUST_NO)) = evList(1, ETIME);
        %如果队列中有顾客
        if (~isempty(qList))
            %队列长度将发生变化，累积前一段时刻的队列长度求和
            nQ = nQ + length(qList)*(evList(1,ETIME) - tQ);
            tQ = evList(1,ETIME);
            %增加繁忙服务台数
            busyN = busyN + 1;
            %添加该顾客的离开事件
            evList(end+1, ETIME) = evList(1, ETIME) + Interval_Serve(qList(1));
            evList(end, ETYPE) = -1;
            evList(end, ECUST_NO) = qList(1);
            %将该顾客移除等待队列
            qList(1) = [];
        end
        %删除该顾客的离开事件
        evList(1,:) = [];
    end
    %对事件列表排序，默认首列升序，即按照时间顺序排序
    evList = sortrows(evList);
end


%各顾客在系统中的逗留时间
t_Wait=t_Leave-t_Arrive;
t_Wait_avg=mean(t_Wait);
%各顾客在系统中的排队时间
t_Queue=t_Wait-Interval_Serve;
t_Queue_avg=mean(t_Queue);

QueLength_avg = nQ/t_Leave(end);

%仿真值与理论值比较
ro = Lambda /(c*Mu);
k=0:(c-1);
P0 = (sum(1./factorial(k).*(Lambda/Mu).^k) +...
    1/factorial(c)/(1-ro)*(Lambda/Mu)^c)^(-1);
Lq = (c*ro)^c*ro/factorial(c)/(1-ro)^2*P0;
Ls = Lq + Lambda/Mu;
Wq = Lq/Lambda;
Ws = Ls/Lambda;

disp(['理论平均逗留时间=',num2str(Ws)]);
disp(['理论平均排队时间=',num2str(Wq)]);
disp(['理论系统中平均队长=',num2str(Lq)]);

disp(['仿真平均逗留时间=',num2str(t_Wait_avg)])
disp(['仿真平均排队时间=',num2str(t_Queue_avg)])
disp(['仿真系统中平均队长=',num2str(QueLength_avg)]);