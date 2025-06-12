%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%first part: compute the arriving time interval, service time
%interval,waiting time, leaving time during the whole service interval
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc
mean_arr = 0.9;
mean_serv = 0.4;
server_num = 3;
cust_num = 100000;
state=zeros(4,cust_num);   
                %state(1,:) represent Arriving time of each customer
                %state(2,:) represent Service time of each customer
                %state(3,:) represent Waiting time of each customer
                %state(4,:) represent Leaving time of each customer 
state(1,:)=exprnd(1/mean_arr,1,cust_num); 
                %arriving time interval between each customer follows exponetial distribution
state(2,:)=exprnd(1/mean_serv,1,cust_num); 
                %service time of each customer follows exponetial distribution
                
for i=1:server_num 
    state(3,1:server_num)=0; 
end 
arr_time=cumsum(state(1,:)); 
                %accumulate arriving time interval to compute arriving time of each customer
state(1,:)=arr_time; 
state(4,1:server_num)=sum(state(:,1:server_num)); 
                %compute living time of first NO.server_num customer by using fomular arriving time + service time
serv_desk=state(4,1:server_num); 
                %create a vector to store leaving time of customers which is in service
for i=(server_num+1):cust_num 
   if arr_time(i)>min(serv_desk) 
       state(3,i)=0;                     
   else  
       state(3,i)=min(serv_desk)-arr_time(i); 
                %when customer NO.i arrives and the server is all busy, the waiting time can be compute by minus arriving time from the minimum leaving time
   end 
   state(4,i)=sum(state(:,i)); 
   for j=1:server_num  
           if serv_desk(j)==min(serv_desk) 
               serv_desk(j)=state(4,i); 
               break 
           end 
                %replace the minimum leaving time by the first waiting customer's leaving time           
   end 
end 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%second part: compute the queue length during the whole service interval
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Timepoint=[state(1,:),state(4,:)];% Arriving time and leaving time
Timepoint=sort(Timepoint);
ArriveFlag=zeros(size(Timepoint));% Flag of arriving time and leaving time
QueLength=zeros(size(Timepoint));% Queue length when each customer arrive or leave
CusNum=zeros(size(Timepoint));% Number of customer when each customer arrive or leav
temp=2;
CusNum(1)=1;
% temp=1;
for i=2:length(Timepoint)
    if (temp<=length(state(1,:)))&&(Timepoint(i)==state(1,temp))
        ArriveFlag(i)=1;
        CusNum(i)=CusNum(i-1)+1;
        temp=temp+1;
    else
        CusNum(i)=CusNum(i-1)-1;
    end
end
for i=1:length(CusNum)
    if CusNum(i)>=server_num+1
        QueLength(i)=CusNum(i)-server_num;
    else
        QueLength(i)=0;
    end
end
Time_interval=zeros(size(Timepoint));% Interval of Events
Time_interval(1)=state(1,1);
for i=2:length(Timepoint)
    Time_interval(i)=Timepoint(i)-Timepoint(i-1);
end
QueLength_avg=sum([0 QueLength].*[Time_interval 0])/sum(Time_interval);% Average length of  queue 

disp(['仿真平均逗留时间=',num2str(mean(state(4,:)-state(1,:)))])
%disp(['仿真平均排队时间=',num2str(mean(state(4,:)-state(2,:)-state(1,:)))])
disp(['仿真平均排队时间=',num2str(mean(state(3,:)))])
disp(['仿真系统中平均队长=',num2str(QueLength_avg)]);
