%Plots the spectra of particles at specified pixels
clear all
%% Parameters
%Specify path to folder with analysis data
WINSD = 'Z:\Zhenyang Jia\03 Experimental\Hyperspectral\250115\150BLZ_750CTR';
sample = 'AuNR_150_750_6analysis';
file = 'AuNR_150_750_6analysis';
cd(strcat(WINSD));
%Specify file name of data you want to plot
load(sample, 'specfin', 'wvlths')
place2save = 'Z:\Zhenyang Jia\03 Experimental\Hyperspectral\250115\150BLZ_750CTR\bareAuNC_5analysis';
mkdir(place2save)

Rawsave = 'AuNR';
Nanosize = '';
colorArea = 'Red'; % Red Orange or Green
Areaval = '';

%Specify pixel coordinates
%automated
% [part_coord,specfin,wvlths] = particle_search(20,WINSD,file,5);
%manual input

part_coord =[
44 13
44 18
36 25
4 32
12 50
5 56
9 57
6 63
25 66
50 58
44 64
51 68
36 77
40 101
6 106
23 110
5 128
26 145
51 145
17 148
42 161
46 175
23 201
40 220
];
temp_coord =part_coord;
part_coord = circshift(part_coord,[0,1]);
matrix_index_ref = part_coord;
matrix_index_ref = circshift(matrix_index_ref,[0,1]);
% specfin = specfin(:,:,141:1000);
specfin_final = sum(specfin(:,:,141:670),3);

%Testing area for shift code
DriftFind='n';
if DriftFind =='y'
    shifter =3;
    for i= 1:length(part_coord)
        y=(part_coord(i,1));
        x=(part_coord(i,2));
        CheckRange = specfin_final(y-shifter:y+shifter,x-shifter:x+shifter);
        maxValue = max(CheckRange(:));
        [y2, x2] = find(CheckRange == maxValue);
        ny = x+(x2-4);
        nx = y+(y2-4);
        part_coord(i,1)=nx;
        part_coord(i,2)=ny;
    end
    newdriftbase=[part_coord(:,2),part_coord(:,1)];
end

high = max(max(specfin_final));
figure1 = figure;
imshow(specfin_final,[0 high])
hold all
for n = 1:size(part_coord,1)
    text(part_coord(n,2),part_coord(n,1),num2str(n),'Color','Green')
end
hgsave('image','-v6')


%Specify cutoff bounds in units of pixels, not nm. 200:1050 gives about 
% 500-850 nm
%NORMAL IS 100 to 670
lower_bound = 1;
upper_bound = 670;
rawwvlths = wvlths(1:670);
wvlths = wvlths(lower_bound:upper_bound);
       
%Describes the number of pixels you want to integrate/particle. Typically 3 or 5.
int_size = 3;                                              
int_var_low = -1*(int_size-1)./2;
int_var_high = (int_size-1)./2;

% n = size(matrix_index_ref,1);
% x_part = matrix_index_ref(n+size(matrix_index_ref,1));
% y_part = matrix_index_ref(n);
% part_spec = zeros(1,1,670);
% for m = int_var_low:int_var_high
%     for l = int_var_low:int_var_high
%         a = x_part + m;
%         b = y_part + l;
%         part_spec = part_spec + specfin(a,b,:);
%     end
% end
% part_spec = squeeze(part_spec);
% background = part_spec(lower_bound:upper_bound);

% % % Remove specfic matrix
% temp = [];
% for i = [1,2,6:15];
%     temp = [temp;matrix_index_ref(i,:)];
% end
% matrix_index_ref = temp;


%% Integrates over the specified square and then plots each specified pixel.
diff = [];
cd(place2save)
all_spec = [];
AllSigyNoi =[];
rawdata = [];
fits =[]; 
for n = 1:size(matrix_index_ref,1)
        % full data background correcting
    x_part = matrix_index_ref(n+size(matrix_index_ref,1));
    y_part = matrix_index_ref(n)+7;
    part_spec = zeros(1,1,670);
    for m = int_var_low:int_var_high
        for l = int_var_low:int_var_high
            a = x_part + m;
            b = y_part + l;
            part_spec = part_spec + specfin(a,b,:);
        end
    end
    part_spec = squeeze(part_spec);
    background = part_spec(1:670);
    
    x_part = matrix_index_ref(n+size(matrix_index_ref,1));
    y_part = matrix_index_ref(n);
    part_spec = zeros(1,1,670);
    for m = int_var_low:int_var_high
        for l = int_var_low:int_var_high
            a = x_part + m;
            b = y_part + l;
            part_spec = part_spec + specfin(a,b,:);
         end
    end
    part_spec = squeeze(part_spec);
    part_specr = part_spec(1:670);
    rawdata = [rawdata,part_specr];
    % fit area
    x_part = matrix_index_ref(n+size(matrix_index_ref,1));
    y_part = matrix_index_ref(n)+7;
    part_spec = zeros(1,1,670);
    for m = int_var_low:int_var_high
        for l = int_var_low:int_var_high
            a = x_part + m;
            b = y_part + l;
            part_spec = part_spec + specfin(a,b,:);
        end
    end
    part_spec = squeeze(part_spec);
    background = part_spec(lower_bound:upper_bound);
    
    x_part = matrix_index_ref(n+size(matrix_index_ref,1));
    y_part = matrix_index_ref(n);
    part_spec = zeros(1,1,670);
    for m = int_var_low:int_var_high
        for l = int_var_low:int_var_high
            a = x_part + m;
            b = y_part + l;
            part_spec = part_spec + specfin(a,b,:);
         end
    end
    part_spec = squeeze(part_spec);
    part_spec = part_spec(lower_bound:upper_bound);
%     part_spec = part_spec - background;
    all_spec = [all_spec,part_spec];
    cd(strcat('Z:\Zhenyang Jia\03 Experimental\Hyperspectral\Hyperspectral Microscope Analysis Code'));
    [param_1,param_2]=fn_lorentz_fit(wvlths',part_spec,1,1);
    a1 = param_1.a1;
    b1 = param_1.b1;
    c1 = param_1.c1;
    resonance = b1;
    FWHM = c1;
    r_list = param_2.rsquare;
    lorentz_fit =(2*a1/pi).*(c1./(4*(wvlths'-b1).^2+c1.^2));
    diff = part_spec-lorentz_fit;
    Noi = std(diff);
    [Notneeded, IndiMax]=min(abs(wvlths-resonance));
    Sigy = part_specr(IndiMax);
    SnN= Sigy/Noi;
    AllSigyNoi = [AllSigyNoi,[Noi;Sigy;SnN]];
    cd(strcat(place2save))
    figure1 = figure;
    hold all
    plot(rawwvlths,part_specr,'b','linewidth',3)
    plot(wvlths,lorentz_fit,'k--','linewidth',3)
    xlabel('Wavelength (nm)','fontsize',32)
    ylabel('Scattering','fontsize',32)
    set(gca,'FontSize',22,'box','on')
    text(0.55,0.9,['\lambda_m_a_x = ',num2str(round(resonance)), ' nm'],'fontsize',20,'Units','normalized')
    text(0.55,0.78,['\Gamma = ',num2str(round(FWHM)), ' nm'],'fontsize',20,'Units','normalized')
    text(0.55,0.66,['S/N = ',num2str(round(SnN))],'fontsize',20,'Units','normalized')
    xlim([500 825])
    ylim([0 max(part_spec)+0.05])
    title(['Particle ',(num2str(n))])
    filename1 = [Nanosize,colorArea,Areaval,num2str(n)];
    saveas(figure1,[filename1,'.tif'])
    hgsave(filename1,'-v6')
    close;
    figure2 = figure;
    hold all
    plot(rawwvlths,part_specr,'b','linewidth',2)
    %plot(wvlths,lorentz_fit,'k--','linewidth',3)
    xlabel('Wavelength (nm)')
    ylabel('Scattering Intensity (arb. u.)')
    title(['Particle ',(num2str(n))])
    set(gca,'box','on')
%    set(gca,'FontSize',22,'box','on')
%    text(0.55,0.9,['\lambda_m_a_x = ',num2str(round(resonance)), ' nm'],'fontsize',20,'Units','normalized')
    %text(0.05,0.78,['\Gamma = ',num2str(round(FWHM)), ' nm'],'fontsize',20,'Units','normalized')
    xlim([400 950])
    ylim([-.02 max(part_spec)+0.02])
    filename2 = [Nanosize,colorArea,Areaval,num2str(n),'r'];
    saveas(figure2,[filename2,'.tif'])
    hgsave(filename2,'-v6')
    temp1 = [wvlths.',part_spec];
    filename3 =[sample,'_','spectra','_',num2str(n),'r'];
%     if Rawsave == 'y'
%         saveas(figure2,[filename2,'.tif'])
%     end
    close;
end
save('spectra','wvlths','all_spec')
