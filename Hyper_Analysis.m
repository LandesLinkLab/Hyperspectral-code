%% Hyper_Analysis
clear all
close all
% tic

%% initial parameters
% cdate = '160823';               % date of data collection
fold = 'Z:\Zhenyang Jia\03 Experimental\Hyperspectral\250115\150BLZ_750CTR\'       %insert search path here

backper = 0.1;                 %fraction of spectra to use in background average
lowercut=50;                   %Pixels to cut from the blue end of spectra
uppercut=50;                   %Pixels to cut from the red side of spectra
lower = 0;                     %lower bound for particle identification (controls how dim of an object is a particle)
upper = 0.5;                    %upper bound for particle identification (controls how bright of an object is selected)
nhood = 1;                      %odd number: nhood by nhood pixels binned to make spectra
rsquarelim = 0.5;
binfac = 'one';                 %binning factor 'one' 'two' or 'four'
bf = 2;                         %1,2,4


%% Opening files, analyzing and saving data in mat  structures

addpath(fold)               
dataloc = fold;
addpath(dataloc)                            % adds directory to path
fid = fopen([dataloc,'\mydata.txt']);       %opens mydata and creates string array
names = textscan(fid, '%s');
fclose(fid);
% names{1,1}(3,1) = cellstr(newname);
nsamp = numel(names{1,1})-2;                %counts number of files to analyze

fcheck=exist([dataloc '\standark.mat'],'file'); %is there a standard already here?
if fcheck==2
    load([dataloc '\standark.mat'])         %if yes, load it
    stanbig = stan;
else
    stan = standardreadindark(strcat(names{1,1}(1),'.tdms'),strcat(names{1,1}(2),'.tdms'),dataloc); %if no, read in the files and make one
end
stan = stan - min(min(min(stan))) + 0.1;    

for c3 = 1:nsamp
    anfunc_lorentz_fit(dataloc,names,c3,stan,backper,lowercut,uppercut,lower,upper,nhood,rsquarelim,binfac)  %analyze all of the files listed in mydata using another function for memory purposes
end
    
%% plotting  -  comment section back in to  see each particle spectrum in the files provided

% for c4=1:nsamp
%     fnew=strcat(names{1,1}(c4+2),'analysis.mat');
%     load(fnew{1,1})
%     srgb=sum(specrgbnorm,3);
% %{
%     [sspec,ind]=sortrows(speccor(:,lowercut+1:end-uppercut),round((1+lowercut)+(numel(wvlths)-uppercut)/2)-lowercut); % 
%     sspec=flipud(sspec);
%     ind=flipud(ind);
%     marksort=mark2(ind,:);
%     sprime=yprime(ind,:);
%     sparams=params(ind,:);
% 
%     pld=figure;
%     set (pld, 'Units', 'normalized', 'Position', [0.1,0.1,0.8,0.8]);
%     subaxis(1,3,1)
%     imshow(specrgbnorm,'InitialMagnification','fit')
%     title('Reconstructed RGB Image')
%     subaxis(1,3,2)
%     imshow(srgb,[0 0.2])
%     title('Grayscale sum image')
%     subaxis(1,3,3)
%     imshow(srgb,[0 0.2])
%     hold on
%     for j = 1:mm
%         text(mark2(j,1),mark2(j,2),num2str(j),'color','r')
%     end
%     hold off
%     title('All Particle Locations')
%   %}  
%     
%  %%   
%         
%      pld=figure;
%         set (pld, 'Units', 'normalized', 'Position', [0.3,0.2,0.25,0.4]);
%     subaxis(2,1,1,'SpacingVertical',0.01,'SpacingHorizontal',0, ...
%             'PaddingLeft',0,'PaddingRight',0,'PaddingTop',0,'PaddingBottom',0, ...
%             'MarginLeft',0.2,'MarginRight',0.05,'MarginTop',0.05,'MarginBottom',0.1,'HoldAxis',1)
%     plot(xbin,yprimesort)
%             set(gca,'fontsize',12,'fontweight','b','Box','on','xticklabel','')
%         axis tight
%             
%     subaxis(2,1,2,'SpacingVertical',0.01,'SpacingHorizontal',0, ...
%             'PaddingLeft',0,'PaddingRight',0,'PaddingTop',0,'PaddingBottom',0, ...
%             'MarginLeft',0.2,'MarginRight',0.05,'MarginTop',0.1,'MarginBottom',0.15,'HoldAxis',1)
%     plot(wvlths(lowercut+1:end-uppercut),medfilt2(speccorsort(:,lowercut+1:end-uppercut),[1 8]))
%     xlabel('wavelength (nm)','fontsize',12,'fontweight','b');
%         set(gca,'fontsize',12,'fontweight','b','Box','on')
%         axis tight
%         y2=ylabel('Intensity (a.u.)','fontsize',12,'fontweight','b');
%         set(y2, 'Units', 'Normalized', 'Position', [-0.12, 1, 0])
% 
%   % histograms of metrics
%         pld=figure;
%         set (pld, 'Units', 'normalized', 'Position', [0.6,0.1,0.2,0.8]);
% 
%    
%         subaxis(3,1,1,'SpacingVertical',0.01,'SpacingHorizontal',0, ...
%             'PaddingLeft',0,'PaddingRight',0,'PaddingTop',0,'PaddingBottom',0, ...
%             'MarginLeft',0.15,'MarginRight',0.1,'MarginTop',0.02,'MarginBottom',0.15,'HoldAxis',1)
%     
%         ecdfplot(paramssort(:,2))
% %        axis([min(paramssort(:,2)) max(paramssort(:,2)) 0 1])
%         text(0.6*(max(paramssort(:,2))-min(paramssort(:,2)))+min(paramssort(:,2)),0.2,['Mean : ', num2str(round(mean(paramssort(:,2))))],'fontsize',12,'fontweight','b');
%         text(0.55*(max(paramssort(:,2))-min(paramssort(:,2)))+min(paramssort(:,2)),0.1,['Std. Dev. : ', num2str(round(std(paramssort(:,2))))],'fontsize',12,'fontweight','b');
%         xlabel('Resonance (nm)','fontsize',12,'fontweight','b');
%         
%         
%         set(gca,'fontsize',12,'fontweight','b','Box','on')
%         
%         y2=ylabel('CDF','fontsize',12,'fontweight','b');
%         set(y2, 'Units', 'Normalized', 'Position', [-0.09, 0.5, 0])
%         
%         subaxis(3,1,2,'SpacingVertical',0.01,'SpacingHorizontal',0, ...
%             'PaddingLeft',0,'PaddingRight',0,'PaddingTop',0,'PaddingBottom',0, ...
%             'MarginLeft',0.15,'MarginRight',0.1,'MarginTop',0.1,'MarginBottom',0.14,'HoldAxis',1)
%         
%         
%         ecdfplot(paramssort(:,3))
% %        axis([min(paramssort(:,3)) max(paramssort(:,3)) 0 1])
% 
%         text(0.6*(max(paramssort(:,3))-min(paramssort(:,3)))+min(paramssort(:,3)),0.2,['Mean : ', num2str(round(mean(paramssort(:,3))))],'fontsize',12,'fontweight','b');
%         text(0.55*(max(paramssort(:,3))-min(paramssort(:,3)))+min(paramssort(:,3)),0.1,['Std. Dev. : ', num2str(round(std(paramssort(:,3))))],'fontsize',12,'fontweight','b');
% 
%         xlabel('FWHM (nm)','fontsize',12,'fontweight','b')
%         ylabel('occurances')
% 
%         set(gca,'fontsize',12,'fontweight','b','Box','on')
%         axis tight
%         y2=ylabel('CDF','fontsize',12,'fontweight','b');
%         set(y2, 'Units', 'Normalized', 'Position', [-0.09, 0.5, 0])
% 
%         subaxis(3,1,3,'SpacingVertical',0.01,'SpacingHorizontal',0, ...
%             'PaddingLeft',0,'PaddingRight',0,'PaddingTop',0,'PaddingBottom',0, ...
%             'MarginLeft',0.15,'MarginRight',0.1,'MarginTop',0.1,'MarginBottom',0.06,'HoldAxis',1)
% 
%         
%         
%         ecdfplot(paramssort(:,1))
%         %axis([min(paramssort(:,1)) max(paramssort(:,1)) 0 1])
%         text(0.6*(max(paramssort(:,1))-min(paramssort(:,1)))+min(paramssort(:,1)),0.2,['Mean : ', num2str(round(mean(paramssort(:,1))))],'fontsize',12,'fontweight','b');
%         text(0.55*(max(paramssort(:,1))-min(paramssort(:,1)))+min(paramssort(:,1)),0.1,['Std. Dev. : ', num2str(round(std(paramssort(:,1))))],'fontsize',12,'fontweight','b');
% 
%         xlabel('Intensity (a.u.)','fontsize',12,'fontweight','b');
%         set(gca,'fontsize',12,'fontweight','b','Box','on')
%         
%         y2=ylabel('CDF','fontsize',12,'fontweight','b');
%         set(y2, 'Units', 'Normalized', 'Position', [-0.09, 0.5, 0])
%         
%         
% % histograms of metrics - horizontal
% %{
%         pld=figure;
%         set (pld, 'Units', 'normalized', 'Position', [0.2,0.3,0.6,0.3]);
% 
%    
%         subaxis(1,3,1,'SpacingVertical',0.01,'SpacingHorizontal',0, ...
%             'PaddingLeft',0,'PaddingRight',0,'PaddingTop',0,'PaddingBottom',0, ...
%             'MarginLeft',0.05,'MarginRight',0.25,'MarginTop',0.1,'MarginBottom',0.2,'HoldAxis',1)
%     
%         ecdfplot(paramssort(:,2))
%         axis([min(paramssort(:,2)) max(paramssort(:,2)) 0 1])
%         text(0.6*(max(paramssort(:,2))-min(paramssort(:,2)))+min(paramssort(:,2)),0.2,['Mean : ', num2str(round(mean(paramssort(:,2))))],'fontsize',12,'fontweight','b');
%         text(0.55*(max(paramssort(:,2))-min(paramssort(:,2)))+min(paramssort(:,2)),0.1,['Std. Dev. : ', num2str(round(std(paramssort(:,2))))],'fontsize',12,'fontweight','b');
%         xlabel('Resonance (nm)','fontsize',12,'fontweight','b');
%         
%         
%         set(gca,'fontsize',12,'fontweight','b','Box','on')
%         
%         y2=ylabel('CDF','fontsize',12,'fontweight','b');
%         set(y2, 'Units', 'Normalized', 'Position', [-0.09, 0.5, 0])
%         
%         subaxis(1,3,2,'SpacingVertical',0.01,'SpacingHorizontal',0, ...
%             'PaddingLeft',0,'PaddingRight',0,'PaddingTop',0,'PaddingBottom',0, ...
%             'MarginLeft',0.15,'MarginRight',0.15,'MarginTop',0.1,'MarginBottom',0.2,'HoldAxis',1)
%         
%         
%         ecdfplot(paramssort(:,3))
%         axis([min(paramssort(:,3)) max(paramssort(:,3)) 0 1])
% 
%         text(0.6*(max(paramssort(:,3))-min(paramssort(:,3)))+min(paramssort(:,3)),0.2,['Mean : ', num2str(round(mean(paramssort(:,3))))],'fontsize',12,'fontweight','b');
%         text(0.55*(max(paramssort(:,3))-min(paramssort(:,3)))+min(paramssort(:,3)),0.1,['Std. Dev. : ', num2str(round(std(paramssort(:,3))))],'fontsize',12,'fontweight','b');
% 
%         xlabel('FWHM (nm)','fontsize',12,'fontweight','b')
% 
%         set(gca,'fontsize',12,'fontweight','b','Box','on')
%         axis tight
%         y2=ylabel('CDF','fontsize',12,'fontweight','b');
%         set(y2, 'Units', 'Normalized', 'Position', [-0.09, 0.5, 0])
% 
%         subaxis(1,3,3,'SpacingVertical',0.01,'SpacingHorizontal',0, ...
%             'PaddingLeft',0,'PaddingRight',0,'PaddingTop',0,'PaddingBottom',0, ...
%             'MarginLeft',0.25,'MarginRight',0.05,'MarginTop',0.1,'MarginBottom',0.2,'HoldAxis',1)
% 
%         
%         
%         ecdfplot(paramssort(:,1))
%         axis([min(paramssort(:,1)) max(paramssort(:,1)) 0 1])
%         text(0.6*(max(paramssort(:,1))-min(paramssort(:,1)))+min(paramssort(:,1)),0.2,['Mean : ', num2str(round(mean(paramssort(:,1))))],'fontsize',12,'fontweight','b');
%         text(0.55*(max(paramssort(:,1))-min(paramssort(:,1)))+min(paramssort(:,1)),0.1,['Std. Dev. : ', num2str(round(std(paramssort(:,1))))],'fontsize',12,'fontweight','b');
% 
%         xlabel('Intensity (a.u.)','fontsize',12,'fontweight','b');
%         set(gca,'fontsize',12,'fontweight','b','Box','on')
%         
%         y2=ylabel('CDF','fontsize',12,'fontweight','b');
%         set(y2, 'Units', 'Normalized', 'Position', [-0.09, 0.5, 0])
%        %} 
% 
%         
%         
%         % R sqaured values
%         plin=figure;
%         
%         ecdfplot(paramssort(:,5))
%         set(gca,'fontsize',10,'fontweight','b','Box','on')
%         %axis([min(paramssort(:,5)) max(paramssort(:,5)) 0 1])
%         
%         plout=figure;
%         set (plout, 'Units', 'normalized', 'Position', [0.3,0.3,0.2,0.3]);
%   
%         ecdfplot(params(:,5))
%         %axis([min(params(:,5)) max(params(:,5)) 0 1])
%         xlabel('Rsquared','fontsize',12,'fontweight','b');
%         set(gca,'fontsize',12,'fontweight','b','Box','on')
%       
%         y2=ylabel('CDF','fontsize',12,'fontweight','b');
%         set(y2, 'Units', 'Normalized', 'Position', [-0.09, 0.5, 0])
% 
%         
%         pli = inset(plout,plin,0.6);
%         text(min(params(:,5)),0.2,['Mean : ', num2str(mean(paramssort(:,5)),2)],'fontsize',12,'fontweight','b');
%         text(min(params(:,5)),0.1,['Std. : ', num2str(std(paramssort(:,5)),2)],'fontsize',12,'fontweight','b');
% 
%         
%        
% %        set (plo, 'Units', 'normalized', 'Position', [0.6,0.1,0.2,0.3]);
% 
%     %%
%     
%     rgbturn(:,:,1)=specrgbnorm(:,:,1)';
%     rgbturn(:,:,2)=specrgbnorm(:,:,2)';
%     rgbturn(:,:,3)=specrgbnorm(:,:,3)';
%     
%       figure; 
%       
%       subaxis(2,1,1,'SpacingVertical',0.01,'SpacingHorizontal',0, ...
%             'PaddingLeft',0,'PaddingRight',0,'PaddingTop',0,'PaddingBottom',0, ...
%             'MarginLeft',0,'MarginRight',0,'MarginTop',0.05,'MarginBottom',0.1,'HoldAxis',1)
%     imshow(rgbturn,'initialmagnification','fit')
%     %title([num2str(newmm),' good, ',num2str(badmm),' bad'])
%     subaxis(2,1,2,'SpacingVertical',0.01,'SpacingHorizontal',0, ...
%             'PaddingLeft',0,'PaddingRight',0,'PaddingTop',0,'PaddingBottom',0, ...
%             'MarginLeft',0,'MarginRight',0,'MarginTop',0.05,'MarginBottom',0.1,'HoldAxis',1)
%     imshow(rgbturn,'initialmagnification','fit')
%     hold on
%     for ki = 1:newmm
%         plot(mark2sort(ki,2),mark2sort(ki,1),'go','MarkerSize',5,'LineWidth',1)
%     end
%     for kj= 1:badmm
%         plot(mark2discard(kj,2),mark2discard(kj,1),'ro','MarkerSize',5,'LineWidth',1)
%     end
%     hold off
%     title([num2str(newmm),' good, ',num2str(badmm),' bad'])
%     
%    
%     for i = 1:1 % 1:newmm %% this can be changed to look at some of the particles per analysis, such as 1:30 to look at the first 30
%         
%         
%         LineWidth=2;
%     FontSize=16;
%         plt=figure;
%         set (plt, 'Units', 'normalized', 'Position', [0.4,0.3,0.3,0.35]);
%         hold on
%         
%         plot(wvlths(1+lowercut:end-uppercut),speccorsort(i,lowercut+1:end-uppercut),'b','LineWidth',LineWidth)
%         plot(wvlths(1+lowercut:end-uppercut),yprimesort(i,:),'r','LineWidth',LineWidth)
%         hold off
%         axis tight
%         box off
%         xlabel('Wavelength (nm)','fontsize',14,'fontweight','b')
%     ylabel('Intensity (a.u.)','fontsize',14,'fontweight','b')
%         set(gca,'fontsize',14,'fontweight','b')
% 
%         %title(['Scattering Spectra and Lorentzian fit, particle ',num2str(i)])
%         axis tight
%         set(gca,'Box','on')
%         labels=cell(4,1);
%         labels{2}=['\lambda_{max} = ' num2str(round(paramssort(i,2))) ' nm'];
%         labels{3}=['FWHM = ' num2str(round(paramssort(i,3))) ' nm'];
%         labels{4}=['R^2 = ' num2str(paramssort(i,5))];
%         text(550,0.7*max(yprimesort(i,:)),labels,'FontSize',14)
%     end
%     rgbturn(:,:,1)=specrgbnorm(:,:,1)';
%     rgbturn(:,:,2)=specrgbnorm(:,:,2)';
%     rgbturn(:,:,3)=specrgbnorm(:,:,3)';
% 
%     figure
%     
%     imshow(sum(rgbturn,3),'InitialMagnification','fit')
%     hold on
%     for j = 1:newmm
%         text(mark2sort(j,2),mark2sort(j,1),num2str(j),'color','g','fontsize',8,'fontweight','b')
%     end
%     hold off
%    
% end
    


%% Plotting for all of the samples
%{
for c4=1:nsamp    
    figure
    subplot(1,2,1)
    imshow(specrgbnorm,'InitialMagnification','fit')
    hold on
    ColOrd = get(gca,'ColorOrder');
    [m,~] = size(ColOrd);
    for ki = 1:mm(1)
        ColRow = rem(ki,m);
        if ColRow == 0
            ColRow = m;
        end
        Col = ColOrd(ColRow,:);
        plot(marker(ki,1),marker(ki,2),'o','MarkerSize',25,'LineWidth',2,'Color',Col)
        hold on
    end
    hold off
    title([names{1,1}(c4+2) 'rgb'])

    subplot(1,2,2)
    plot(wvlths(1,lowercut+1:end),spfl{c4}(:,lowercut+1:end))
    xlabel('wavelength (nm)'); ylabel('Scattering Intensity (a.u.)')
    hold on
    plot(wvlths(1,lowercut+1:end),yprime2{c4}(:,:),':')
    hold off
    title('local bg and flatfield corrected')
 
end
%% 

for c6 = 1:mm(1)
    figure
    subplot(1,2,1)
    imshow(specrgbnorm{c4},'InitialMagnification','fit')
    hold on
    plot(marker(c6,1),marker(c6,2),'o','MarkerSize',25,'LineWidth',2)
    hold off
    subplot(1,2,2)
    hold on
    for c7= 1:nsamp
        ColRow = rem(c7,m);
        if ColRow == 0
            ColRow = m;
        end
        Col = ColOrd(ColRow,:);
        plot(wvlths(1,lowercut+1:end),spfl{c7}(c6,lowercut+1:end),'Color',Col)
        text(700,max(max(spfl{1}(c6,lowercut+1:end)))-c7*(((max(max(spfl{1}(c6,lowercut+1:end)))-min(spfl{1}(c6,lowercut+1:end))))/nsamp),names{1,1}(c7+2),'Color', Col)
    end    
    hold off
    xlabel('wavelength (nm)'); ylabel('Scattering Intensity (a.u.)')
end

%%
figure
hold on
for c5 = 1:nsamp
    ColRow = rem(c5,m);
    if ColRow == 0
        ColRow = m;
    end
    Col = ColOrd(ColRow,:);
    plot(wvlths(1,lowercut+1:end),spfl{c5}(:,lowercut+1:end),'Color',Col)
    text(500,6-c5*.45,names{1,1}(c5+2),'Color', Col)
end
hold off
title('local bg and flatfield corrected - data')
xlabel('wavelength (nm)'); ylabel('Scattering Intensity (a.u.)')

figure
hold on
for c5 = 1:nsamp
    ColRow = rem(c5,m);
    if ColRow == 0
        ColRow = m;
    end
    Col = ColOrd(ColRow,:);
    plot(wvlths(1,lowercut+1:end),yprime2{c5}(:,:),'Color',Col)
    text(500,6-c5*.45,names{1,1}(c5+2),'Color', Col)
end
hold off
title('local bg and flatfield corrected - fits')
xlabel('wavelength (nm)'); ylabel('Scattering Intensity (a.u.)')
    %}
%}
% totalt=toc