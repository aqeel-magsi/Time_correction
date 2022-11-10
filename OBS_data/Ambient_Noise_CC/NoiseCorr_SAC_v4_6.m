%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extracting cross-correlation function (CF) from cross correlation of daily 
% ambient seismic noise data. 
% The current version deals with daily long sac format data and output daily 
% cross-correlation function (in non-overlapping period bands) in mat format 
% and the final stacked broadband CFs in ascii format The input data must 
% have the same sampling frequency.

% In this version, we do one-bit or temporal normalization cross-correlation 
% for different bands seperately. And then normalize the daily CFs in different 
% bands and stack them together to form the broadband CFs. This may help to 
% improve SNR of CFs in different bands than the normal one broadband prcoessing.
% Before using this code, please read the main  function NoiseCorrMBand_SAC_v4 
% very carefully. You have to change seisfile1 and seisfile2 in the main 
% function to read data successfully.

% - by Huajian Yao, 2014 Dec 23, USTC
% hjyao@ustc.edu.cn   huajianyao@gmail.com

% Reference: 
% Yao, H., van der Hilst R.D., and de Hoop, M.V., 2006. Surface-wave array 
% tomography in SE Tibet from ambient seismic noise and two-station analysis 
% : I - Phase velocity maps. Geophys.J. Int., Vol. 166(2), 732-744, 
% doi: 10.1111/j.1365-246X.2006.03028.x. 
% Yao, H., Gouedard, P., McGuire, J., Collins, J. and van der Hilst, R.D., 
% 2011. Structure of young East Pacific Rise lithosphere from ambient noise 
% correlation analysis of fundamental- and higher-mode Scholte-Rayleigh waves, 
% Comptes Rendues Geoscience de l'Acad��mie des Sciences, 343��571?583��
% doi:10.1016/j.crte.2011.04.004


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function NoiseCorr_SAC_v4

%%---------Define the value of essential parameters in the programs---------
%__________________________________________________________________________
StaFile = 'obs.loc';   % file containing the station information
RespFile = 'OBSstaResp.txt'; % file containing corresponding station response file (RESP. format)
IndexRmResp = 0; % index of performing instrument response removal or not: =1 remove, will use RespFile; otherwise, RespFile is not used
datadir = '/home/aqeel/Research_PhD/OBS_Data/Noise_Correlation/data_hyd_6h_org/' ; % folder containing all sac data for all stations (data must be named with Julian day ...)
PeriodBand = [2 5]; % array of period bands of the CFs, unit: sec; you can define multiple non-overlapping bands here,
%PeriodBand = [2 10];                         % e.g.: PeriodBand = [5 10; 10 20; 20 40], and you will get CFs in three period bands,and also in the broad band 5-40s
%fsNew = 10; % data resampling frequency; original sample freq / fsNew must be an integer!!!
fsNew = 10;

IndexWhiteSpec = 1; % = 1, spectrum whitening; otherwise use original spectrum
indexCorrMethod = 2; % = 1: one-bit cross correlation; = 2: temporal normalization cross-correlation [suggest to use this one]
yearrange = [2018 2019]; % year range of the data
dayrange = [1 366]; % range of days in a year for CF calculation 
hourrange = [1 4]; % range of days in a year for CF calculation 
%hourrange = [1 2]; % range of hours in a day for CF calculation 
MaxLagTime = 200;  % maximum time length (in seconds) for each side of CFs
SegHour = 2; % data segment length (in hour) for data preprocessing (removing response, spectral whitenning)
cmp1 = 'hyd'; % data component (Z, E, N) for the first station for the processing
cmp2 = 'hyd'; % data component (Z, E, N) for the second station for the processing
outCFdir = '/home/aqeel/Research_PhD/OBS_Data/Noise_Correlation/data_hyd_6h_org/NCCF_org/'; % output folder for CFs
%__________________________________________________________________________
mkdir(outCFdir); % create output CF folder

%% read data file
allsta = struct('name', {}, 'lat', {}, 'lon', {}, 'elev',{}, 'respfile', {});
i=0;
fstat = fopen(StaFile,'r');
while ~feof(fstat)
    name = fscanf(fstat,'%s',1);
    if ~strcmp(name,'')
        i=i+1;
	   allsta(i).name = name;  %station name
%	   allsta(i).net = fscanf(fstat,'%s',1) %station network
	   allsta(i).lat = fscanf(fstat,'%f',1); %station latitude
	   allsta(i).lon = fscanf(fstat,'%f',1); %station longitude      
       allsta(i).elev = fscanf(fstat,'%f',1); %station elevation
    else
        break
    end
    temp = fgetl(fstat);
end
NumSta = i;

%% read instrument file if we want to remove instrument response
if IndexRmResp == 1
    fresp = fopen(RespFile,'r');
    for i = 1:NumSta
        name = fscanf(fresp,'%s',1);
        allsta(i).respfile = fscanf(fresp,'%s',1);
    end
end


%% frequency band info
NumPB = size(PeriodBand,1); % number of period band
Tmax = max(PeriodBand(:,2)); % max period for filtering
Tmin = min(PeriodBand(:,1)); % min period for filtering
TRange = [0.75*Tmin Tmin Tmax 1.5*Tmax]; % period range for initial band pass filtering
freqrange = 1./TRange(end:-1:1); % freq range for initial band pass filtering , reverse the Trange: 
                       % [f_low_cut f_low_pass f_high_pass f_high_cut]
                       
%% create folders to store mat and ascii format CFs

matCFdir = [outCFdir cmp1 '-' cmp2 '/mat']; % output data folder for daily CFs in different bands
asciiCFdir = [outCFdir cmp1 '-' cmp2 '/ascii']; % output data folder for stacked broadband CFs
[s,mess,messid] = mkdir(matCFdir); % create station pair folder for saving mat CF data files
[s,mess,messid] = mkdir(asciiCFdir); % create station pair folder for saving ascii CF data files

%% cross-correlation for every combination of two stations (also include
%% autocorrelation functions
    
%% loop for first station
for i = 1:NumSta-1

    %% loop for second station
	for j = i+1 : NumSta  
        
        display(['station pair =  ' allsta(i).name '    ' allsta(j).name]);

        stainfo(1) = allsta(i);
        stainfo(2) = allsta(j);
        MaxShiftNum = round(MaxLagTime*fsNew);
        CFtime = ((-MaxShiftNum):1:MaxShiftNum)'/fsNew;
        nptCF = length(CFtime);
        CFdata = struct('year',0, 'day',0, 'NCF', zeros(nptCF, NumPB));
        
        
        ncf = 0; 
        tic
        %% loop for year
        for year= yearrange(1):yearrange(2);
            yearstr = num2str(year);
            %% loop for day
            for day = dayrange(1):dayrange(2) % loop for each julian day in a year

                % obtain the string for the day
                if day < 10
                    daystr = ['00' num2str(day)]
                elseif day >= 10 && day < 100
                    daystr = ['0' num2str(day)]
                elseif day >= 100
                    daystr = num2str(day);
                end
	         %% loop for hour
                for hour = hourrange(1):hourrange(2) % loop for each julian day in a year
			    hourstr = num2str(hour);
                
                % specify the data file name (NOTE: USER HAS TO MODIFY THIS PART TO HAVE CORRECT INPUT DATA FILE !!!!!!!)
                % -------------------------------------------------------------------------------------------------------
                seisfile1 = [datadir allsta(i).name '/' yearstr '/' daystr  '_' hourstr  '_' cmp1 '.sac'];
                seisfile2 = [datadir allsta(j).name '/' yearstr '/' daystr  '_' hourstr  '_' cmp2 '.sac'];
                %seisfile1 = [datadir allsta(i).name '/' yearstr '/' daystr '_' cmp1 '.sac'];
               % seisfile2 = [datadir allsta(j).name '/' yearstr '/' daystr '_' cmp2 '.sac'];
                % -------------------------------------------------------------------------------------------------------
                
                
                 display(['read data ...' yearstr '  '  daystr]);
                % tic
                % read the specified sac data file, if this sac file does not exit,
                % the sac structure is empty (i.e., isempty(sta1) = 1)
                sta1 = readsac(seisfile1); % use this if readsacFS is not working
                sta2 = readsac(seisfile2); % use this if readsacFS is not working
                % toc
                
                % sta1 = readsacFS(seisfile1,0); 
                % sta2 = readsacFS(seisfile2,0);
                 
                % check whether the daily data is existing and longer than one hour
                sta1_dataok = (~isempty(sta1) && sta1.NPTS>3600 *1/sta1.DELTA);
                sta2_dataok = (~isempty(sta2) && sta2.NPTS>3600 *1/sta2.DELTA);

        
                if  sta1_dataok && sta2_dataok % sac data file for station 1 & 2 exists
                    
                    display(['Now processing: ' yearstr '  '  daystr '   ...... ']);  
            
                    if sta1.DELTA ~= sta2.DELTA
                        display('Error: data sampling rate is not the same!');
                        break;
                    end
                    
                    % resampling sta1 waveform to a new sampling frequency
                    if (1/sta1.DELTA) > fsNew
                        DecimateR_org = (1/sta1.DELTA)/fsNew;
                        DecimateR = round(DecimateR_org);
                        if abs(DecimateR_org - DecimateR) > 0.001
                            display('Error: resampling frequency!');
                            break;
                        end
                        nn = floor(sta1.NPTS/DecimateR);
                        if (nn*DecimateR+1) <= sta1.NPTS
                            ReSampleWave = decimate(sta1.DATA1(1:nn*DecimateR+1), DecimateR);
                        else
                            ReSampleWave = decimate([sta1.DATA1(1:nn*DecimateR); sta1.DATA1(nn*DecimateR)], DecimateR);
                        end   
                        sta1.DELTA = 1.0/fsNew;
                        sta1.NPTS = length(ReSampleWave);
                        sta1.DATA1 = ReSampleWave;
                        clear ReSampleWave
                    end   
                    
                    % resampling sta2 waveform to a new sampling frequency
                    if (1/sta2.DELTA) > fsNew
                        DecimateR_org = (1/sta2.DELTA)/fsNew;
                        DecimateR = round(DecimateR_org); 
                        if abs(DecimateR_org - DecimateR) > 0.001
                            display('Error: resampling frequency!');
                            break;
                        end
                        nn = floor(sta2.NPTS/DecimateR);
                        if (nn*DecimateR+1) <= sta2.NPTS
                            ReSampleWave = decimate(sta2.DATA1(1:nn*DecimateR+1), DecimateR);
                        else
                            ReSampleWave = decimate([sta2.DATA1(1:nn*DecimateR); sta2.DATA1(nn*DecimateR)], DecimateR);
                        end    
                        sta2.DELTA = 1.0/fsNew;
                        sta2.NPTS = length(ReSampleWave);
                        sta2.DATA1 = ReSampleWave;
                        clear ReSampleWave
                    end                     
                    
                    % display('preprocess data:');
                    % tic
                    % create band pass filter
                    SampleF = 1/sta1.DELTA;
                    LowF = (2/SampleF)*freqrange(1);
                    HighF = (2/SampleF)*freqrange(end);
                    [B, A] = butter(2,[LowF, HighF]);
                    
                    sta1.DATA1 = detrend(sta1.DATA1 - mean(sta1.DATA1)); % demean and detrend the data
                    % figure(1); subplot(2,1,1); hold off; plot(sta1.DATA1(1:1000)); 
                    sta1.DATA1 = filtfilt(B, A, sta1.DATA1); % initially bandpass filter the data
                    % figure(1); subplot(2,1,2); hold off; plot((1:1000)/SampleF,sta1.DATA1(1:1000)); waitforbuttonpress

                    % bandpass filter the waveform and remove instrument response
                    % (if IndexRmResp ~= 1, do not remove instrument response)
                    bfwave1 = rmResp_bpfilter(sta1.DATA1, 1/sta1.DELTA, freqrange, allsta(i).respfile, IndexRmResp, IndexWhiteSpec,  SegHour);  
                    sta1.DATA1 = bfwave1; 
                    clear bfwave1;
                 %   display(['   station 1 *** ' allsta(i).name ' ***']);
                    
                    sta2.DATA1 = detrend(sta2.DATA1 - mean(sta2.DATA1)); % detrend the data
                    % figure(1); subplot(2,1,1); hold off; plot(sta2.DATA1(1:1000));
                    sta2.DATA1 = filtfilt(B, A, sta2.DATA1); % initially bandpass filter the data
                    % figure(1); subplot(2,1,2); hold off; plot((1:1000)/SampleF,sta2.DATA1(1:1000)); waitforbuttonpress
                   
                    % function call : read and remove instrument response; also band-pass filtering the wave trains
                    bfwave2 = rmResp_bpfilter(sta2.DATA1, 1/sta2.DELTA, freqrange, allsta(j).respfile, IndexRmResp, IndexWhiteSpec, SegHour); 
                    sta2.DATA1 = bfwave2; 
                    clear bfwave2;
                 %   display(['         station 2 *** ' allsta(j).name ' ***']);
      
                    stadist = deg2km(distance(allsta(i).lat, allsta(i).lon, allsta(j).lat, allsta(j).lon)) % station distance (km)
                    % toc
                    
                    % display('noise cross-correlation:');
                    % tic
                    % function call : multiband noise cross-correlation with time-domain normalization 
                    CFcnMB = CrossCorrelation(sta1, sta2, PeriodBand, MaxLagTime, indexCorrMethod);
                    %[CFcnMB,lags]= CrossCorrelation(sta1, sta2, PeriodBand, MaxLagTime, indexCorrMethod);
                    if sum(sum(isnan(CFcnMB))) == 0 % no NaN data in CFcn
                        ncf = ncf + 1;
                        CFdata(ncf).NCF = CFcnMB;
                        CFdata(ncf).year = year; 
                        CFdata(ncf).day = day;
                       % CFdata(ncf).hour = hour;
			
                    end
                    % toc
                    
                end % end if: data existing
               end %for hour
            end % end for loop: day
         
            
        end % % end for loop: year
        
        %% stack daily CF (normalized) of each freq band, then stack each band to form a broadband CF
        %  save CFdata & CF info to mat files; save stacked broadband CF
        %  data to ascii files which can be used for later dispersion
        %  analysis
        
        if ncf > 0
            cross_info = [allsta(i).lon,allsta(i).lat, allsta(i).elev; allsta(j).lon, allsta(j).lat, allsta(j).elev];
            stackCF = zeros(length(CFtime),NumPB+1);
            nnday = 0;
            for nd = 1:ncf
                dailyCFmax = max(CFdata(nd).NCF);
                if max(dailyCFmax) > 0
                    stackCF(:,1:NumPB) = stackCF(:,1:NumPB) + CFdata(nd).NCF*diag(1./max(CFdata(nd).NCF)); % stack normalized daily CF for each period band
                    nnday = nnday + 1;
                end
            end
            stackCF(:,1+NumPB) = sum(stackCF(:,1:NumPB), 2); % sum each period band CF to obtain the broadband CF
            stackCF = stackCF*diag(1./max(stackCF)); % normalized CF in each frequency band
            
           CFcnmatfile = [matCFdir '/' cmp1 cmp2 '_' allsta(i).name '-' allsta(j).name '.mat']; % output matCF data file name
           CFcnascfile = [asciiCFdir '/' cmp1 cmp2 '_' allsta(i).name '-' allsta(j).name '_' num2str(nnday) 'd.dat']; % ascii file name for saving broadband CF
    
            
            save(CFcnmatfile, 'stainfo', 'yearrange', 'dayrange', 'PeriodBand', 'SampleF', 'MaxLagTime', 'cmp1', 'cmp2');    % save sta and noise CF info
            save(CFcnmatfile, 'CFtime', 'stackCF', 'nnday', '-append');
            save(CFcnmatfile, 'CFdata',  '-append');
            
            CFcn = zeros(MaxShiftNum+1, 3);
            CFcn(:,1) = CFtime((MaxShiftNum+1):(2*MaxShiftNum+1));
            CFcn(:,2) = stackCF((MaxShiftNum+1):(2*MaxShiftNum+1),1+NumPB);
            CFcn(:,3) = stackCF((MaxShiftNum+1):-1:1,1+NumPB);
            save(CFcnascfile,'cross_info','CFcn','-ASCII');  % save station info and CF (broadband stack) to an ascii file
            
        end
        toc
        clear CFdata CFtime stackCF stainfo cross_info CFcn
        
    end % end for loop: the first station    
end % end for loop: the second station

clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---                Read Instrument Response File                         ---%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Normfactor, Numzeros, Numpoles, Respzero, Resppole] = Rd_InstruRespFile(RespFile)
%%
fname = fopen(RespFile,'r');
% skip 1 - 18 line
for i = 1:18
temp1 = fgetl(fname);
end

%read line 19
temp1 = fscanf(fname,'%s',4);
Normfactor = fscanf(fname,'%f',1);
temp1 = fgetl(fname);

%read line 20
temp1 = fgetl(fname);
%read line 21
temp1 = fscanf(fname,'%s',4);
Numzeros = fscanf(fname,'%f',1);
temp1 = fgetl(fname);
%read line 22
temp1 = fscanf(fname,'%s',4);
Numpoles = fscanf(fname,'%f',1);
temp1 = fgetl(fname);
%read line 23, 24: zeroes header
temp1 = fgetl(fname);
temp1 = fgetl(fname);
%read zeros
Respzer o = zeros(1, Numzeros);
for i = 1:Numzeros
   temp1 = fscanf(fname,'%s',1);
   temp = fscanf(fname,'%d',1);
   realpart = fscanf(fname,'%e',1);
   imagpart = fscanf(fname,'%e',1);
   Respzero(i) = complex(realpart, imagpart);
   temp1 = fgetl(fname);
end
%read 2 lines: poles header
temp1 = fgetl(fname);
temp1 = fgetl(fname);
%read poles
Resppole = zeros(1, Numpoles);
for i = 1:Numpoles
   temp1 = fscanf(fname,'%s',1);
   temp = fscanf(fname,'%d',1);
   realpart = fscanf(fname,'%e',1);
   imagpart = fscanf(fname,'%e',1);
   Resppole(i) = complex(realpart, imagpart);
   temp1 = fgetl(fname);
end
fclose(fname);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% remove instrument responses, whiten, and band-pass filter the wave 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data Processing: remove instrument response, band-pass filtering
function bfwave = rmResp_bpfilter(seisdata, fs, freqrange, respfile, IndexRmResp, IndexWhiteSpec, SegHour)
% seisdata: data vector
% fs: sampling frequency
% freqrange = [f_low_cut f_low_pass f_high_pass f_high_cut] % filter freq band
%              f_low_cut < f_low_pass < f_high_pass < f_high_cut, e.g.,
%              [0.005 0.01 5 10] Hz
% respfile: response file name
% IndexRmResp = 1 : read RESP.* response file and remove response
%             else: do not remove instrument response
% IndexWhiteSpec: = 1, spectrum whitening; otherwise keep the original spectrum
% SegHour: data segment length (in hour) for processing the data
%          (remove instrument response, whitening, filtering, 

% define instrument response information
resp = struct('Amp',0,...
    'NumZeros',0,...
    'NumPoles',0,...
    'Zeros',0,...
    'Poles',0,...
    'poly_num',0,...
    'poly_den',0);

%% 
HighF = freqrange(3);  % high pass freq
LowF = freqrange(2);   % low pass freq
%% HighF must be less than fs/2, otherwise set to be 0.99*(fs/2)
if HighF > (fs/2)
    HighF = 0.99*(fs/2);
    display(['HighF error! HighF is reset to ', num2str(HighF)]);
end

%% read instrument response files
if IndexRmResp == 1 % read RESP.* response file
    [resp.Amp, resp.Numzeros, resp.Numpoles, resp.Zeros, resp.Poles] = Rd_InstruRespFile(respfile);
    [resp.poly_num, resp.poly_den] = zp2tf(resp.Zeros', resp.Poles', resp.Amp);
end

LowFMin = max(freqrange(1), 0); % lowest freq for response removal
HighFMax = min(freqrange(4), fs/2); % highest freq for reponse removal

SegLength = round(SegHour*3600*fs); % length of data parts for fft and processing (e.g., whitening)
if mod(SegLength,2) == 1 % ensure the even number of SegLength
    SegLength = SegLength + 1;
end
staNumPt = length(seisdata); % number of points in data
fftnum = ceil(staNumPt/SegLength); % number of segments for the data

%%
bfwave=zeros(1,staNumPt);

for k=0:(fftnum-1)
    %%
    if k ~=(fftnum-1)
        datapart = detrend(seisdata((1 + k*SegLength):(k + 1)*SegLength));
        fftlength=SegLength;
        fftdata=fft(datapart,fftlength);
        clear datapart
    else
        datapart = detrend(seisdata((1 + k*SegLength):staNumPt));
        fftlength = 2^(nextpow2(staNumPt - k*SegLength));
        fftdata=fft(datapart, fftlength);
        clear datapart
    end

    if max(abs(fftdata)) > 0 && fftlength > fs*3000 % to prevent the processing of zero or NaN data in the input data segment, more than half hour  data
        fftdata = reshape(fftdata, 1, fftlength);
        %%
        f(1:(fftlength/2+1)) = fs*(0:(fftlength/2))/fftlength;
        delta_f = fs/fftlength;


        MinFPoint = max(2, ceil(LowFMin/delta_f));
        MaxFPoint = min(fftlength/2, floor(HighFMax/delta_f));

        % remove instrument response
        if IndexRmResp == 1
            % remove instrument response: the first half frequency spectrum
            w(1:(fftlength/2+1)) = 2*pi*f(1:(fftlength/2+1));
            h = freqs(resp.poly_num, resp.poly_den, w); % obtain instrument response
            h = h/max(abs(h)); % normalize the instrument response
            nn =  MinFPoint:MaxFPoint;
            % Y = XH -> X = Y/H -> X = Y*conj(H)/abs(H)^2
            h = h/max(abs(h)); % normalize the amplitude of the intrument response
            fftdata(nn) = fftdata(nn).*conj(h(nn))./(abs(h(nn)).^2 + 0.01);  % water-level deconvolution
            fftdata(1:MinFPoint) = 0;
            fftdata(MaxFPoint:(fftlength/2+1)) = 0;
            fftdata((fftlength/2+2):fftlength)=conj(fftdata((fftlength/2):-1:2)); % treat another half spectrum
        end

        %% spectrum whitenning (divide the smooth spectrum amplitude from running average)
        if IndexWhiteSpec == 1
            nn =  MinFPoint:MaxFPoint;
            datafamp = abs(fftdata(nn));
            winsize = max(round(0.02/delta_f), 11);
            if mod(winsize,2) == 0
                winsize = winsize + 1;
            end
            shiftpt = round((winsize+1)/2);
            datafampnew = [ones(1,winsize)*datafamp(1) datafamp ones(1,winsize)*datafamp(end)];
            datafampsmooth = filter(ones(1,winsize)/winsize,1,datafampnew);
            datafampsmooth2 = datafampsmooth((winsize+shiftpt):(winsize+shiftpt+length(nn)-1));
            KK = find(datafampsmooth2 > 0);
            JJ = isnan(datafampsmooth2);
            fftdata(nn(JJ)) = 0;
            fftdata(nn(KK)) = fftdata(nn(KK))./datafampsmooth2(KK) ; 
            
    %         figure(1); hold off; subplot(2,1,1); plot(datafamp); hold on; plot(datafampsmooth2, 'r');
    %         subplot(2,1,2); plot(abs(fftdata(nn)));
        end

        %% band pass filtering
        LowPtN = round(LowF/delta_f);
        HighPtN = round(HighF/delta_f); 
        nptdfs = round((LowF - LowFMin)/delta_f);
        if nptdfs >= 4
            nn = (LowPtN - nptdfs):(LowPtN-1);
            taperwin = hann(2*nptdfs-1)';
            fftdata(1:(LowPtN - nptdfs -1))=0; 
            % figure(99); hold off; subplot(2,1,1); hold off; plot(abs(fftdata(nn)),'r');
            fftdata(nn) = taperwin(1:nptdfs).*fftdata(nn);
            % hold on; plot(abs(fftdata(nn)),'b--');
        end

        nptdfs = round((HighFMax - HighF)/delta_f);
        nn = (HighPtN + 1):(HighPtN + nptdfs);
        if nptdfs >= 4
            taperwin = hann(2*nptdfs-1)';
            % subplot(2,1,2); hold off; plot(abs(fftdata(nn)),'r');
            fftdata(nn)= taperwin(nptdfs:end).*fftdata(nn);   
            fftdata((HighPtN + nptdfs + 1):(fftlength/2+1)) = 0;
            % hold on; plot(abs(fftdata(nn)),'b--');
        end

        fftdata((fftlength/2+2):fftlength)=conj(fftdata((fftlength/2):-1:2));

        %% time domain filtered data
        band_data_T=real(ifft(fftdata));

    %     figure(2); hold off; subplot(2,1,1), semilogx(f(1:(fftlength/2+1)), abs(fftdata(1:(fftlength/2+1)))); xlim([LowFMin HighFMax]);
    %     subplot(2,1,2), plot(band_data_T);
    %     waitforbuttonpress

        if k~=(fftnum-1)
            bfwave((1 + k*SegLength):(k + 1)*SegLength)=band_data_T;
        else
            bfwave((1 + k*SegLength):staNumPt)=band_data_T(1:(staNumPt - k*SegLength));
        end
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% multibands cross-correlation function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CFcnMB = CrossCorrelation(sta1, sta2, PeriodBand, MaxLagTime, indexCorrMethod)
% sta1 & sta2: input sac data structure
% PeriodBand: period bands [StartT1 EndT1; StartT2 EndT2; ...] for cross-correlation 
%             (first filter the data in each band, then normalization, then cross correlation)
% MaxLagTime: cross correlation time lag: between [-MaxLagTime MaxLagTime] s
% indexCorrMethod: cross-correlation method: =1 for one-bit ; =2 for temporal normalization 
% CFcnMB: output CF matrix for multibands

MonthDays = [31 29 31 30 31 30 31 31 30 31 30 31; 31 28 31 30 31 30 31 31 30 31 30 31];

% one bit cross correlation
SampleF = round(1/sta1.DELTA);
SampleT = sta1.DELTA;

if SampleF == 0
    error('Error! Sampling Frequency of data should not be zero!');
end

% if the first point of the two wain trains start at differnet year (e.g., 2003, 2004)
% Here we only deal with the case of one year difference
if (sta2.NZYEAR - sta1.NZYEAR) == 1
    if mod(sta1.NZYEAR, 4) == 0  
        YearDays = 366;
    else
        YearDays = 365;  
    end            
    sta2.NZJDAY = sta2.NZJDAY + YearDays;
elseif (sta1.NZYEAR - sta2.NZYEAR) == 1  
    if mod(sta2.NZYEAR, 4) == 0  
        YearDays = 366;
    else
        YearDays = 365;  
    end            
    sta1.NZJDAY = sta1.NZJDAY + YearDays;
elseif abs(sta1.NZYEAR - sta2.NZYEAR) > 1
    display('Year Error! ');
    exit
end

DeltaTInitial = (sta2.NZJDAY - sta1.NZJDAY)*24*3600 + (sta2.NZHOUR - sta1.NZHOUR)*3600 + (sta2.NZMIN - sta1.NZMIN)*60 + ... 
                (sta2.NZSEC - sta1.NZSEC) + (sta2.NZMSEC - sta1.NZSEC)/1000; 
DeltaTInitial = round(DeltaTInitial*SampleF);

% let the first point of the two wave trains start at same time (the later one)
if DeltaTInitial > 0
    if (sta1.NPTS-DeltaTInitial) > 0
        sta1.DATA1(1:(sta1.NPTS-DeltaTInitial)) = sta1.DATA1((DeltaTInitial+1):sta1.NPTS);
        sta1.DATA1((sta1.NPTS-DeltaTInitial+1):sta1.NPTS) = 0;    
    else
        sta1.DATA1(1:sta1.NPTS) = 0;
    end
    DeltaTInitial = 0;
elseif DeltaTInitial < 0
    DeltaTInitial = abs(DeltaTInitial);
    if (sta2.NPTS-DeltaTInitial) > 0
        sta2.DATA1(1:(sta2.NPTS-DeltaTInitial)) = sta2.DATA1((DeltaTInitial+1):sta2.NPTS);
        sta2.DATA1((sta2.NPTS-DeltaTInitial+1):sta2.NPTS) = 0; 
    else
        sta2.DATA1(1:sta2.NPTS) = 0;
    end
    DeltaTInitial = 0;
end

PointNum = min(sta2.NPTS, sta1.NPTS);
MaxShiftNum = round(MaxLagTime/SampleT);
MinShiftNum = -MaxShiftNum;
ShiftNum = MaxShiftNum - MinShiftNum + 1;  % also ShiftNum = 2*MaxTravT*SampleF + 1

NumPB = size(PeriodBand,1); % number of period band

% GFcnMB = zeros(NumPB, MaxShiftNum+1, 3); % multiband Green's function: output
% CFcnMB = zeros(NumPB, MaxShiftNum+1, 3); % multiband cross-correlation function: output

CFcnMB = zeros(ShiftNum, NumPB);

%% noise cross-correlation for multi-frequency bands
for np = 1:NumPB

    % obtain the Start and End period for the period band
    StartT = PeriodBand(np, 1);
    EndT = PeriodBand(np, 2);
    % bandpass filter and one-bit normalize the waveform 
    LowF = (2/SampleF)/EndT;
    HighF = (2/SampleF)/StartT;
    [B, A] = butter(2,[LowF, HighF]);
    seisdata1 = filtfilt(B, A, sta1.DATA1);
    seisdata2 = filtfilt(B, A, sta2.DATA1);
    if indexCorrMethod == 1  % one-bit normalization
        seisdata1 = sign(seisdata1);
        seisdata2 = sign(seisdata2);
    elseif indexCorrMethod == 2 % temporal normalization by dividing the smooth amplitude (from running average)
        winsize = round(EndT*2*SampleF);
        if mod(winsize,2) == 0
            winsize = winsize + 1;
        end
        shiftpt = round((winsize+1)/2);

        tempamp = [ones(1,winsize)*abs(seisdata1(1)) abs(seisdata1) ones(1,winsize)*abs(seisdata1(end))];
        tempamp = filter(ones(1,winsize)/winsize, 1, tempamp); % running average of amplitude
        tempamp2 = tempamp((shiftpt+winsize):(shiftpt+winsize+sta1.NPTS-1));
        KK = find(tempamp2 >0);
        JJ = isnan(tempamp2);
         seisdata1(JJ) = 0;
         %figure; hold off; subplot(2,1,1); plot(seisdata1(1:10000)); hold on; plot(tempamp((shiftpt+winsize):(shiftpt+winsize+10000-1)),'r');
        seisdata1(KK) = seisdata1(KK)./tempamp2(KK); %where the NaN will occur?
        %subplot(2,1,2); plot(seisdata1(1:10000));
        clear tempamp tempamp2 KK JJ

        tempamp = [ones(1,winsize)*abs(seisdata2(1)) abs(seisdata2) ones(1,winsize)*abs(seisdata2(end))];
        tempamp = filter(ones(1,winsize)/winsize, 1, tempamp); % running average of amplitude
        tempamp2 = tempamp((shiftpt+winsize):(shiftpt+winsize+sta2.NPTS-1));
        KK = find(tempamp2 >0);
        JJ = isnan(tempamp2);
        seisdata2(JJ) = 0;
         %figure; hold off; subplot(2,1,1); plot(seisdata2(1:10000)); hold on; plot(tempamp((shiftpt+winsize):(shiftpt+winsize+10000-1)),'r');
        seisdata2(KK) = seisdata2(KK)./tempamp2(KK);
        % subplot(2,1,2); plot(seisdata2(1:10000));
        clear tempamp tempamp2 KK JJ
    end

     %display('Onebit normalization:');  
    % one-bit cross-correlation
    OneBitCross = xcorr(seisdata2(1:PointNum), seisdata1(1:PointNum), MaxShiftNum);

     %display('Onebit Cross correlation:');

    % band-pass filter the correlation function in [LowF, HighF]
    OneBitCross = filtfilt(B, A, OneBitCross);

%     % Differentiate the correlation function to get the Green's fcn
%     GreenFcn = zeros(ShiftNum,1);
% 
%     for i = 2:(ShiftNum-1)
%         GreenFcn(i) = (OneBitCross(i+1) - OneBitCross(i-1))/2;
%     end
%     GreenFcn(1) = OneBitCross(2) - OneBitCross(1);
%     GreenFcn(ShiftNum) = OneBitCross(ShiftNum) - OneBitCross(ShiftNum - 1);
%     % if max(GreenFcn>0)
%     %     GreenFcn = GreenFcn/max(GreenFcn);
%     % end
%     
%     GFcn = zeros(MaxShiftNum+1, 3); % Green's function    
%     GFcn(:, 1) = (0:1:MaxShiftNum)'/SampleF;
%     GFcn(:, 2) = - GreenFcn((MaxShiftNum+1):ShiftNum);
%     GFcn(:, 3) = GreenFcn((MaxShiftNum+1):-1:1);
% 
%     CFcn = zeros(MaxShiftNum+1, 3); % cross correlation function
%     CFcn(:, 1) = (0:1:MaxShiftNum)'/SampleF;
%     CFcn(:, 2) = OneBitCross((MaxShiftNum+1):ShiftNum);
%     CFcn(:, 3) = OneBitCross((MaxShiftNum+1):-1:1);
    
%     GFcnMB(np,:,:) = GFcn;
%     CFcnMB(np,:,:) = CFcn;
       
    CFcnMB(:,np) = OneBitCross;
    CFtime = ((-MaxShiftNum):1:MaxShiftNum)'/SampleF;
    
end
