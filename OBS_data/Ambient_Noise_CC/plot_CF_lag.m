% plot the CFs.
% ZHU Gaohua, June 2018, CUHK
dmat = dir('*.mat'); % only looking for .mat Files
N = length(dmat);    % number of .mat Files
bpn=1; % PeriodBand = [2 5; 5 10];

%%
for ii = 1
    load(['hydhyd_' num2str(ii) '.mat'])
    time = CFtime(2100:2350)';
    %time = CFtime(1:end)';
    firstday = 1;
    lengday = length(CFdata);
    lastday=lengday;
    bd=[2018 09 07];
    ed=[2019 10 15];
    sg=12;
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
        for j = -2:2 %%
            k=i+j;
            if k > lengday; k = lengday; end
            if k <1; k = 1; end
            zz = CFdata(k).NCF(2100:2350,1)'./max(CFdata(k).NCF(2100:2350,1))/(abs(j)+1);
            %zz = CFdata(k).NCF(1:end,1)'./max(CFdata(k).NCF(1:end,1))/(abs(j)+1);
            zz2= zz2+zz;

        end

        z2(i,:)= zz2;
        z2s=z2s+z2(i,:);
    end
   
     z2s=z2s/max(abs(z2s));
     fig=figure(ii);
     h1=subplot(3,1,1);
     imagesc(time,day,z2);
     colormap(jet)
     %convert hour segments into date, bd here is begin day
     YTickStr = char(datetime(bd, 'InputFormat', 'yyyy/MM/dd', 'Format', 'dd/MM/yy') + hours(day*sg));
     %yval = ylim(gca);
     L=round(linspace(1,size(YTickStr,1),5));
     set(gca,'YTick',day(:,L),'YTickLabel', YTickStr(L,:));
     xlim([-100 100]);
     
    h3=subplot(3,1,2);
    %set(h3,'position',[0.1 0.32 0.3 0.1])
    plot(time,z2s,'k','linewidth',1)
    xlabel('Lag time (s)');
    xlim([-100 100]);
    %set(fig,'PaperOrientation','landscape');
    
        
%for ii = 1:3
    
%load(['bhzbhz_' num2str(ii) '.mat'])
a1=[];
b1=[];
%time = CFtime(1250:1500)';

for i=1:lengday
    [ref_value,index11]=max(z2(8,:)); %get first peak value and index as reference 
    [max_value,index]=max(z2(i,:)); %max peak value and index for all days
    d=time(index); %extract time values w.r.t to index
    d11=time(index11); %extract time for first reference
    df_value=minus(max_value,ref_value); %substract and get the difference between peaks
    d_ind= minus(d,d11); %substract and get the difference between lag time w.r.t to first day 
    %d=[x,y]
    a1(i,:)=[i,d_ind];
end
h4=subplot(3,1,3);
    d=1:lengday;  %total no of days
    % Get coefficients of a line fit through the data.
    coefficients = polyfit(a1(:,1),a1(:,2),1);
    % Create a new x axis with exactly 1000 points (or whatever you want).
    xFit = linspace(min(d), max(d), lengday);
    X=xFit';
    % Get the estimated yFit value for each of those 1000 new x locations.
    yFit = polyval(coefficients ,xFit);
    Y=yFit';
    b1=[X,Y];
    % Plot everything.
    plot(a1(:,1),a1(:,2), 'b.', 'MarkerSize', 15); % Plot training data.
    hold on; % Set hold on so the next plot does not blow away the one we just drew.
    plot(b1(:,1),b1(:,2), 'r-', 'LineWidth', 2); % Plot fitted line.
    grid on
    hold off
    L=round(linspace(1,size(YTickStr,1),5));
    set(gca,'XTick',day(:,L),'XTickLabel', YTickStr(L,:));
    
    %xlabel('No. of Segments (6 hours)');
    ylabel('Lag time(s)');
    xlim([1 lengday]);
    ylim([-25 10]);
    %set(fig,'PaperOrientation','landscape');
    F = sprintf('lag_data_%d.txt',ii);
    dlmwrite(F, a1, 'delimiter',' ')
    K = sprintf('linear_fit_%d.txt',ii);
    dlmwrite(K, b1, 'delimiter',' ')
    S = sprintf('Date_%d.txt',ii);
    dlmwrite(S,YTickStr,'delimiter',' ')
    %hold off
%dlmwrite('lag_data'num2str(ii) '.txt',a1,'delimiter',' ')
print(fig,'-dpdf','-fillpage');
 clear  CFtime z2s a1 b1 X Y xFit yFit day
end




