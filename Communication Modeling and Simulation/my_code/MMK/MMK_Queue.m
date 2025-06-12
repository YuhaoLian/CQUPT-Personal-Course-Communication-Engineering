function [t_Wait_avg,t_Queue_avg,QueLength_avg] = MMK_Queue(Quene_Para)

Sim_Customer_Num = Quene_Para.Sim_Customer_Num;
Lambda=Quene_Para.Lambda;
Mu=Quene_Para.Mu;
K = Quene_Para.K;
Time_Arrive=Quene_Para.Time_Arrive;
Time_Leave=Quene_Para.Time_Leave;
Interval_Arrive=Quene_Para.Interval_Arrive;
Interval_Serve=Quene_Para.Interval_Serve;

%按照分布，初始化到达间隔和服务间隔
Interval_Arrive=-log(rand(1,Sim_Customer_Num))/Lambda;
Interval_Serve=-log(rand(1,Sim_Customer_Num))/Mu;
%通过累计求和，计算顾客的达到时刻；
Time_Arrive = cumsum(Interval_Arrive);

%定义时间列表的三个参量的索引值：
%事件时刻，事件类型（1代表到达，-1代表离开），客户编号
ETIME = 1; ETYPE = 2; ECUST_NO = 3;
%初始化事件列表，加入第一个客户到达的事件
evList(1,ETIME) = Time_Arrive(1);
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
        if (busyN < K)
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
        if (nextN <= Sim_Customer_Num)
            evList(end+1, ETIME) = Time_Arrive(nextN);
            evList(end, ETYPE) = 1;
            evList(end, ECUST_NO) = nextN;
            nextN = nextN + 1;
        end
        %离开事件
    else
        %繁忙服务台减一
        busyN = busyN - 1;
        %记录该顾客的离开事件
        Time_Leave(evList(1, ECUST_NO)) = evList(1, ETIME);
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
t_Wait=Time_Leave-Time_Arrive;
t_Wait_avg=mean(t_Wait);
%各顾客在系统中的排队时间
t_Queue=t_Wait-Interval_Serve;
t_Queue_avg=mean(t_Queue);

QueLength_avg = nQ/Time_Leave(end);

%仿真值与理论值比较
ro = Lambda /(K*Mu);
k=0:(K-1);
P0 = (sum(1./factorial(k).*(Lambda/Mu).^k) +...
    1/factorial(K)/(1-ro)*(Lambda/Mu)^K)^(-1);
Lq = (K*ro)^K*ro/factorial(K)/(1-ro)^2*P0;
Ls = Lq + Lambda/Mu;
Wq = Lq/Lambda;
Ws = Ls/Lambda;

figure;
subplot(1,2,1);
plot(1:Sim_Customer_Num,t_Wait);
xlabel('顾客数(number)')
ylabel('逗留时间(t/s)')
title('顾客在系统中的逗留时间')
subplot(1,2,2);
plot(1:Sim_Customer_Num,t_Queue);
xlabel('顾客数(number)')
ylabel('排队时间(t/s)')
title('各顾客在系统中的排队时间')
ylim([0 max(t_Queue)+2]);

disp(['理论平均逗留时间=',num2str(Ws)]);
disp(['理论平均排队时间=',num2str(Wq)]);
disp(['理论系统中平均队长=',num2str(Lq)]);

disp(['仿真平均逗留时间=',num2str(t_Wait_avg)])
disp(['仿真平均排队时间=',num2str(t_Queue_avg)])
disp(['仿真系统中平均队长=',num2str(QueLength_avg)]);
end

