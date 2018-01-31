% CODED BY : PUNEET DHEER
% DATE : 11-01-2018
% INPUT:
% data = multichannel each column is a channel (input is column matrix)
% Ws = window size in sample point
% Lp = Shifting the window by Lp sample point
% freq_band = [0.1 4; 4 7; 7 13; 13 30; 30 80] indvidual band spectral entropy
% sfreq = sampling frequency (in sample points)
% entire_band = [0.1 4 7 13 30 80] grouped together
% fmax = for realtive entropy computation 1 or 0 for not
% fmaxx = 1 for in upto frequency band or 0 for upto sfreq/2
%
% OUTPUT:
% Entire_RANGE_spentropy = total spectral entropy [1 last_band]
% feat_cell = within band entropy
% Entire_spentropy = entire band group wise with relative power SPECTRAL ENTROPY
% Eigenchannels = correlation between channels
% Eigenfreqbands = correlation between bands

function [Eigenchannels, Eigenfreqbands , Entire_RANGE_spentropy, feat_cell, Entire_spentropy]=spectral_entropy_v1(data, Ws, Lp, freq_band, sfreq, entire_band, fmax,fmaxx)

Lw=1;
Z=Ws;

windows=ceil((length(data)-Ws+1)/Lp); 

% entire_band=horzcat(freq_band(1,1),freq_band(:,2)');

No_freq_bands=size(freq_band,1);
feat_cell=cell(1,size(data,2)); % for each channel

for i = 1 : size(data,2)
    feat_cell{1,i}=NaN(windows,No_freq_bands); %freq band for each channel cell window wise
end


for i=1:windows
    
    fprintf('T_Windows: %d \n',i);
    
    window_data=data(Lw:Z,:); %windowed data
    window_data = zscore(window_data,[],1);
    LEN=length(window_data);
    % PSD = (AMPLITUDE SPECTRUM).^2 
    fftx = power(abs(fft(window_data)/(LEN/2)),2); %extract the magnitude by using fft or pwelch

%     fftx = (abs(fft(window_data)).^2)*(1/(sfreq*LEN));
%     [fftx,f] = pwelch(window_data,2048,0,2048,sfreq);  N_freq=2048;
        
%     fftx = realsqrt(fftx); % AMPLITUDE SPECTRUM
   
    fftx(1,:) = 0; %remove DC component which is 0 Hz in FFT
    
    N = length(fftx);  % length of psd
    
    N_freq=N;
    
    mid = floor(N/2);
        
    fftx = fftx(1:mid+1,:); % only positive frequency (one sided frequency)
    
%     fftx = fftx ./ ( N );
    
    N = length(fftx);
    
    
    scale=(N_freq/sfreq);
    
%     frequency_Hz=(0:(N-1))./scale;
    
%     plot(frequency_Hz,fftx)
    
    
%%%%%%%%%%%%%%%%%%%% SPECTRAL ENTROPY [1 last_band]%%%%%%%%%%%%%%%%%

if fmaxx==1 
     fftxx = fftx(ceil(freq_band(1))+1:freq_band(end)+1,:);
 end


fftx_ENTIRE = bsxfun(@rdivide,fftxx,sum(fftxx)); %normalize the fftx in range [0 1] for PDF
   
       
Entire_RANGE_spentropy(i,:) = -sum(fftx_ENTIRE.*log(fftx_ENTIRE+eps))/ log(size(fftx_ENTIRE,1));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





    
%%%%%%%%%%%%%%%%%%%-- within-band SPECTRAL ENTROPY--%%%%%%%%%%%

    for band=1:No_freq_bands
        
        if(band==1)
            start = ceil(freq_band(band,1)*scale);
        else
            start = band_range(end)-1;
        end
        
        band_range = start : floor(freq_band(band,2)*scale);
        band_range = band_range+1;
        band_range(band_range<1) = 1; 
        band_range(band_range>N) = N;    
        
        pr=  bsxfun(@rdivide,fftx(band_range,:),sum(fftx(band_range,:)));  %[0 1]
        
        SpecEn= -sum( pr.*log(pr) )./log(size(pr,1)); %Normalized [0 1]
        
        for ii = 1 : size(data,2)
             feat_cell{1,ii}(i,band) = SpecEn(ii) ;
        end   
       
    end   
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    




%%%%%%%%%%-- entire band group wise with relative power SPECTRAL ENTROPY--%%%%%%

 freq_seg_band = round(size( window_data,1)/sfreq*entire_band)+1; %frequency window 

 if fmax==1 % for the relative in the max frequency range
%      fftx = fftx(1:80+1,:);
     fftx = fftx(freq_seg_band(1):freq_seg_band(end),:);  
 end

SCALE=(freq_seg_band(1)-1);
  
fftx = bsxfun(@rdivide,fftx,sum(fftx)); %normalize the fftx in range [0 1] for PDF

% fftx = fftx(freq_seg_band(1):freq_seg_band(end),:);  

spect = zeros(length(entire_band)-1,size(window_data,2)); %frequency based spectrum
    
    for n=1:length(entire_band)-1
        
        spect(n,:) = sum(fftx(freq_seg_band(n)-SCALE:freq_seg_band(n+1)-SCALE,:));
        
    end
    
Entire_spentropy(i,:) = -sum(spect.*log(spect))/ log(size(spect,1));


    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Eigenchannels(i) = max(abs(CORRELATION(spect))); %between channels
Eigenfreqbands(i) = max(abs(CORRELATION(spect.')));  %between bands

    Lw=Lw+Lp;
    Z=Z+Lp;
        
    
end




end
