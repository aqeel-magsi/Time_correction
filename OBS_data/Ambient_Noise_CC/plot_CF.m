% plot the CFs.
% ZHU Gaohua, June 2018, CUHK
% Aqeel plot and calculate  the lag time, May 2021 
dmat = dir('*.mat'); % only looking for .mat Files
N = length(dmat);    % number of .mat Files
bpn=1; % PeriodBand = [2 5; 5 10];

%%
for ii = 1:N
    load(['bhzbhz_' num2str(ii) '.mat'])
    time = CFtime(1:end)';
    firstday = CFdata(1).day;
    lengday = length(CFdata);
    lastday=CFdata(lengday).day;
    bd=[2020 04 18];
    ed=[2020 06 12];
    d1 = datenum(bd);d2 = datenum(ed);
    day = firstday:firstday+lengday-1;
    %day1=firstday:366;
    %day2=1:lastday;
    %day=[day1,day2];
    [x,y]=meshgrid(time,day);
    z1 = zeros(size(x));
    z2 = zeros(size(x));
    z3 = zeros(size(x));
    z1s = zeros(size(time));
    z2s = zeros(size(time));
    z3s = zeros(size(time));
    
    for i = 1:lengday
        zz2 = 0;
        zz=0;
        for j = 0:0  %%
            k=i+j;
            if k > lengday; k = lengday; end
            if k <1; k = 1; end
            zz = CFdata(k).NCF(1:end,1)'./max(CFdata(k).NCF(1:end,1))/(abs(j)+1);
            %zz = sta(k).data(1:end,bpn)'./max(sta(k).data(1:end,bpn))/(abs(j)+1);
            zz2= zz2+zz;

        end

        z2(i,:)= zz2;
        z2s=z2s+z2(i,:);
    end
   
    z2s=z2s/max(abs(z2s));
    fig=figure(ii);
    h1=subplot(3,1,1);
    %set(h1,'position',[0.1 0.48 0.3 0.45])
    imagesc(time,day,z2)
    colormap(jet)
    xlim([-10 10]);
    
    NumTicks = 4;
    L = get(gca,'YLim');
    set(gca,'YTick',linspace(L(1),L(2),NumTicks));
    datetick('y','mmmdd','keeplimits','keepticks');
    ylabel('(Date)');
    
%     ax=gca;
%     ax.XGrid = 'on';
    
    %figure(4);plot(time,z2s)
    h3=subplot(3,1,2);
    %set(h3,'position',[0.1 0.32 0.3 0.1])
    plot(time,z2s,'k','linewidth',1)
    xlabel('Lag time (s)');
    xlim([-10 10]);
    %set(fig,'PaperOrientation','landscape');
    
    clear  CFdata CFtime z2s
%% Calculate the lag time 
for ii = 1:3
    
load(['bhzbhz_' num2str(ii) '.mat'])
a1=[] 
time = CFtime(1:end)';

for i=1:lengday
    [ref_value,index1]=max(z2(1,:)); %get first peak value and index as reference 
    [max_value,index]=max(z2(i,:)); %max peak value and index for all days
    %[x,y]=find(z2(i,:)==maximum); 
    %lags_correspond_tomaximum_value
    d=time(index); #extract time values w.r.t to index
    d11=time(index11); #extract time for first reference
    df_value=minus(max_value,ref_value); %substract and get the difference between peaks
    d_ind= minus(d11,d); %substract and get the difference between lag time w.r.t to first day 
    %d=[x,y]
    a1(i,:)=[d_ind,df_value];
end
h4=subplot(3,1,3);
    x=1:55  %total no of days
    plot(x,a1(:,1),'k','linewidth',0.6) %plot days and lag times
    xlabel('No. of Days');
    ylabel('Lag time(s)');
    xlim([1 60]);
    ylim([-2 3]);
    %set(fig,'PaperOrientation','landscape');
dlmwrite('lag_data.txt',a1,'delimiter',' ')
print(fig,'-dpdf','-fillpage');
end
end


