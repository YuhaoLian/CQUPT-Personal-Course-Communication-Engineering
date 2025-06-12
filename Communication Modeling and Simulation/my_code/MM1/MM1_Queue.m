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

Interval_Arrive=-log(rand(1,Sim_Customer_Num))/Lambda;%����ʱ����
Interval_Serve=-log(rand(1,Sim_Customer_Num))/Mu;%����ʱ��
Time_Arrive(1)=Interval_Arrive(1);%�˿͵���ʱ��

ArriveNum(1)=1;
for i=2:Sim_Customer_Num
    Time_Arrive(i)=Time_Arrive(i-1)+Interval_Arrive(i);
    ArriveNum(i)=i;
end

Time_Leave(1)=Time_Arrive(1)+Interval_Serve(1);%�˿��뿪ʱ��
LeaveNum(1)=1;
for i=2:Sim_Customer_Num
    if Time_Leave(i-1)<Time_Arrive(i)
        Time_Leave(i)=Time_Arrive(i)+Interval_Serve(i);
    else
        Time_Leave(i)=Time_Leave(i-1)+Interval_Serve(i);
    end
    LeaveNum(i)=i;
end

t_Wait=Time_Leave-Time_Arrive; %���˿���ϵͳ�еĶ���ʱ��
t_Wait_avg=mean(t_Wait);
t_Queue=t_Wait-Interval_Serve;%���˿���ϵͳ�е��Ŷ�ʱ��
t_Queue_avg=mean(t_Queue);

Timepoint=[Time_Arrive,Time_Leave];%ϵͳ���뿪�¼��������¼�������ʱ��
Timepoint=sort(Timepoint);
ArriveFlag=zeros(size(Timepoint));%����ʱ���־
CusNum=zeros(size(Timepoint));%��������뿪�¼�����ʱϵͳ�Ĺ˿�����
temp=2;
CusNum(1)=1;
for i=2:length(Timepoint)
    if (temp<=length(Time_Arrive))&&(Timepoint(i)==Time_Arrive(temp))%�����������κ��¼�������ͬʱ�����ļ���
        CusNum(i)=CusNum(i-1)+1;
        temp=temp+1;
        ArriveFlag(i)=1;
    else
        CusNum(i)=CusNum(i-1)-1;
    end
end
%�����κ������¼���ʱ����
Time_interval=zeros(size(Timepoint));
Time_interval(1)=Time_Arrive(1);
for i=2:length(Timepoint)
    Time_interval(i)=Timepoint(i)-Timepoint(i-1);
end
%ϵͳ��ƽ���˿�������
CusNum_fromStart=[0 CusNum];
CusNum_avg=sum(CusNum_fromStart.*[Time_interval 0] )/Timepoint(end);

%�����¼�����ʱϵͳ���еĳ���
QueLength=zeros(size(CusNum));
for i=1:length(CusNum)
    if CusNum(i)>=2
        QueLength(i)=CusNum(i)-1;%ϵͳ��ֻ��1������̨
    else
        QueLength(i)=0;
    end
end
QueLength_avg=sum([0 QueLength].*[Time_interval 0] )/Timepoint(end);%ϵͳƽ���ӳ�

%����ֵ������ֵ�Ƚ�
% disp(['����ƽ������ʱ��=',num2str(1/(Mu-Lambda))]);
% disp(['����ƽ���Ŷ�ʱ��=',num2str(Lambda/(Mu*(Mu-Lambda)))]);
% disp(['����ϵͳ��ƽ���ӳ�=',num2str(Lambda*Lambda/(Mu*(Mu-Lambda)))]);
% 
% disp(['����ƽ������ʱ��=',num2str(t_Wait_avg)])
% disp(['����ƽ���Ŷ�ʱ��=',num2str(t_Queue_avg)])
% disp(['����ϵͳ��ƽ���ӳ�=',num2str(QueLength_avg)]);
end

