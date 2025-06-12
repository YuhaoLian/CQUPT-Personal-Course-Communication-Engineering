clear;clc;close all;
% M/M/K�Ŷ�ϵͳ����
% �����¼��б�ģ��

SimTotal=100000; %����˿�����
Lambda=0.9; %������Lambda
Mu=0.4; %������Mu��
c = 3; %����̨����

%���շֲ�����ʼ���������ͷ�����
Interval_Arrive=-log(rand(1,SimTotal))/Lambda;
Interval_Serve=-log(rand(1,SimTotal))/Mu;
%ͨ���ۼ���ͣ�����˿͵Ĵﵽʱ�̣�
t_Arrive = cumsum(Interval_Arrive);

%��ʼ���뿪ʱ��
t_Leave=zeros(1,SimTotal);
%����ʱ���б����������������ֵ��
%�¼�ʱ�̣��¼����ͣ�1�����-1�����뿪�����ͻ����
ETIME = 1; ETYPE = 2; ECUST_NO = 3;
%��ʼ���¼��б������һ���ͻ�������¼�
evList(1,ETIME) = t_Arrive(1);
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
        if (busyN < c)
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
        if (nextN <= SimTotal)
            evList(end+1, ETIME) = t_Arrive(nextN);
            evList(end, ETYPE) = 1;
            evList(end, ECUST_NO) = nextN;
            nextN = nextN + 1;
        end
    %�뿪�¼�
    else
        %��æ����̨��һ
        busyN = busyN - 1;
        %��¼�ù˿͵��뿪�¼�
        t_Leave(evList(1, ECUST_NO)) = evList(1, ETIME);
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
t_Wait=t_Leave-t_Arrive;
t_Wait_avg=mean(t_Wait);
%���˿���ϵͳ�е��Ŷ�ʱ��
t_Queue=t_Wait-Interval_Serve;
t_Queue_avg=mean(t_Queue);

QueLength_avg = nQ/t_Leave(end);

%����ֵ������ֵ�Ƚ�
ro = Lambda /(c*Mu);
k=0:(c-1);
P0 = (sum(1./factorial(k).*(Lambda/Mu).^k) +...
    1/factorial(c)/(1-ro)*(Lambda/Mu)^c)^(-1);
Lq = (c*ro)^c*ro/factorial(c)/(1-ro)^2*P0;
Ls = Lq + Lambda/Mu;
Wq = Lq/Lambda;
Ws = Ls/Lambda;

disp(['����ƽ������ʱ��=',num2str(Ws)]);
disp(['����ƽ���Ŷ�ʱ��=',num2str(Wq)]);
disp(['����ϵͳ��ƽ���ӳ�=',num2str(Lq)]);

disp(['����ƽ������ʱ��=',num2str(t_Wait_avg)])
disp(['����ƽ���Ŷ�ʱ��=',num2str(t_Queue_avg)])
disp(['����ϵͳ��ƽ���ӳ�=',num2str(QueLength_avg)]);