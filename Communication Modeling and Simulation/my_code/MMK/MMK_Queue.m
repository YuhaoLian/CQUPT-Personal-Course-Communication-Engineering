function [t_Wait_avg,t_Queue_avg,QueLength_avg] = MMK_Queue(Quene_Para)

Sim_Customer_Num = Quene_Para.Sim_Customer_Num;
Lambda=Quene_Para.Lambda;
Mu=Quene_Para.Mu;
K = Quene_Para.K;
Time_Arrive=Quene_Para.Time_Arrive;
Time_Leave=Quene_Para.Time_Leave;
Interval_Arrive=Quene_Para.Interval_Arrive;
Interval_Serve=Quene_Para.Interval_Serve;

%���շֲ�����ʼ���������ͷ�����
Interval_Arrive=-log(rand(1,Sim_Customer_Num))/Lambda;
Interval_Serve=-log(rand(1,Sim_Customer_Num))/Mu;
%ͨ���ۼ���ͣ�����˿͵Ĵﵽʱ�̣�
Time_Arrive = cumsum(Interval_Arrive);

%����ʱ���б����������������ֵ��
%�¼�ʱ�̣��¼����ͣ�1�����-1�����뿪�����ͻ����
ETIME = 1; ETYPE = 2; ECUST_NO = 3;
%��ʼ���¼��б������һ���ͻ�������¼�
evList(1,ETIME) = Time_Arrive(1);
evList(1,ETYPE) = 1;
evList(1,ECUST_NO) = 1;
%��ʼ������
qList = [];
%��һ������Ŀͻ���ţ�
nextN = 2;
%��æ����̨����
busyN = 0;
%����ʱ���Ȩ��͵Ķ��г���
nQ = 0;
%���г�����һ�α仯��ʱ��
tQ = 0;
%���߼������¼��б�ǿ���Ϊѭ������
while(~isempty(evList))
    %�����¼�
    if (evList(1,ETYPE) == 1)
        %����з���̨����
        if (busyN < K)
            %��æ����̨��һ
            busyN = busyN + 1;
            %���øÿͻ����뿪�¼������߼�Ӧ����ɾ�������¼���������뿪�¼�
            %�˴��Ż�Ϊֱ���޸��¼����ͺ�ʱ�̣�
            evList(1,ETIME ) = evList(1,ETIME) + Interval_Serve(evList(1,ECUST_NO));
            evList(1,ETYPE) = -1;
            %����̨û�п���
        else
            %���г��Ƚ������仯���ۻ�ǰһ��ʱ�̵Ķ��г������
            nQ = nQ + length(qList)*(evList(1,ETIME) - tQ);
            tQ = evList(1,ETIME);
            %���˿��������
            qList(end+1) = evList(1, ECUST_NO);
            %ɾ���ù˿ʹﵽ�¼�
            evList(1,:) = [];
        end
        %������пͻ���Ҫ��������¹˿͵����¼��������¼��б�
        if (nextN <= Sim_Customer_Num)
            evList(end+1, ETIME) = Time_Arrive(nextN);
            evList(end, ETYPE) = 1;
            evList(end, ECUST_NO) = nextN;
            nextN = nextN + 1;
        end
        %�뿪�¼�
    else
        %��æ����̨��һ
        busyN = busyN - 1;
        %��¼�ù˿͵��뿪�¼�
        Time_Leave(evList(1, ECUST_NO)) = evList(1, ETIME);
        %����������й˿�
        if (~isempty(qList))
            %���г��Ƚ������仯���ۻ�ǰһ��ʱ�̵Ķ��г������
            nQ = nQ + length(qList)*(evList(1,ETIME) - tQ);
            tQ = evList(1,ETIME);
            %���ӷ�æ����̨��
            busyN = busyN + 1;
            %��Ӹù˿͵��뿪�¼�
            evList(end+1, ETIME) = evList(1, ETIME) + Interval_Serve(qList(1));
            evList(end, ETYPE) = -1;
            evList(end, ECUST_NO) = qList(1);
            %���ù˿��Ƴ��ȴ�����
            qList(1) = [];
        end
        %ɾ���ù˿͵��뿪�¼�
        evList(1,:) = [];
    end
    %���¼��б�����Ĭ���������򣬼�����ʱ��˳������
    evList = sortrows(evList);
end


%���˿���ϵͳ�еĶ���ʱ��
t_Wait=Time_Leave-Time_Arrive;
t_Wait_avg=mean(t_Wait);
%���˿���ϵͳ�е��Ŷ�ʱ��
t_Queue=t_Wait-Interval_Serve;
t_Queue_avg=mean(t_Queue);

QueLength_avg = nQ/Time_Leave(end);

%����ֵ������ֵ�Ƚ�
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
xlabel('�˿���(number)')
ylabel('����ʱ��(t/s)')
title('�˿���ϵͳ�еĶ���ʱ��')
subplot(1,2,2);
plot(1:Sim_Customer_Num,t_Queue);
xlabel('�˿���(number)')
ylabel('�Ŷ�ʱ��(t/s)')
title('���˿���ϵͳ�е��Ŷ�ʱ��')
ylim([0 max(t_Queue)+2]);

disp(['����ƽ������ʱ��=',num2str(Ws)]);
disp(['����ƽ���Ŷ�ʱ��=',num2str(Wq)]);
disp(['����ϵͳ��ƽ���ӳ�=',num2str(Lq)]);

disp(['����ƽ������ʱ��=',num2str(t_Wait_avg)])
disp(['����ƽ���Ŷ�ʱ��=',num2str(t_Queue_avg)])
disp(['����ϵͳ��ƽ���ӳ�=',num2str(QueLength_avg)]);
end

