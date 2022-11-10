% plot the CFs.
% ZHU Gaohua, June 2018, CUHK

clear all
close all
clc

dmat = dir('*.mat'); % only looking for .mat Files
N = length(dmat);    % number of .mat Files
bpn=1; % PeriodBand = [2 5; 5 10];

%%
for ii = 1:N
    load(['hydhyd_' num2str(ii) '.mat'])
    time = CFtime(1:end)';
    firstday = CFdata(1).day;
    lengday = length(CFdata);
    lastday=CFdata(lengday).day;
    bd=[2018 10 01];
    ed=[2019 10 10];
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
        for j = -5:5   %%
            k=i+j;
            if k > lengday; k = lengday; end
            if k <1; k = 1; end
            %zz = CFdata(k).NCF(1:end,1)'./max(CFdata(k).NCF(1:end,1))/(abs(j)+1);
            zz = CFdata(k).NCF(1:end,bpn)'./max(CFdata(k).NCF(1:end,bpn));
            zz2= zz2+zz;

        end

        z2(i,:)= zz2;
        z2s=z2s+z2(i,:);
    end
   
    z2s=z2s/max(abs(z2s));
    fig=figure(ii);
    h1=subplot(2,1,1);
    %set(h1,'position',[0.1 0.48 0.3 0.45])
    imagesc(time,day,z2)
    colormap(jet)
    xlim([-200 200]);
    
    NumTicks = 6;
    L = get(gca,'YLim');
    set(gca,'YTick',linspace(L(1),L(2),NumTicks));
    datetick('y','yymmm','keeplimits','keepticks');
    ylabel('Date');
    
%     ax=gca;
%     ax.XGrid = 'on';
    
    %figure(4);plot(time,z2s)
    h3=subplot(2,1,2);
    %set(h3,'position',[0.1 0.32 0.3 0.1])
    plot(time,z2s,'k','linewidth',1)
    xlabel('Lag time (s)');
    xlim([-200 200]);
    %set(fig,'PaperOrientation','landscape');
    print(fig,'-dpdf','-fillpage');
    
    clear  CFdata CFtime z2s
    
end


