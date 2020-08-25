clear all
close all
clc

%% Input information
% Check the correction factors and field resolution

filenum = input('File name? (e.g. iAM382  the file extension is implicitly .csv)\n', 's');
Dye = input('Type of dyes? (e.g. D for Donor, A for Acceptor and B for Both)\n', 's');

color = input('Which dyes? (C for Cy3-Cy5, A for a488-a594 or D for deconvolution of a488-a594)\n', 's');

if strcmpi(color, 'C')
    color={'Cy5'; 'tmr'; 'tmrfret'};
elseif strcmpi(color, 'A')
    color={'a594'; 'a488'; 'a594fret'};
elseif strcmpi(color, 'D')
    color={'dw_a594'; 'dw_a488'; 'dw_a594fret'};
end

Correction = input('Correct images? (y/n)\n', 's');

if strcmpi(Correction, 'y')
    CorrDonor = input('Correction ratio for Donor intensity? (e.g. for Cy3 is 0.278, a488 is 0.085 and dw_a488 is 0.0605)\n');
    CorrAcceptor = input('Correction ratio for Acceptor intensity? (e.g. for Cy5 is 0.09, a594 is 0.056 and dw_594 is 0.0637)\n');
end

disp('If there is a group of dots at a distance lower than 2 µm, only the brigher dot is considered.')
disp('Dots with aberrant intensity values are removed. (e.g. 10000 when normal is 10)')
MinDots = input('Minimum dots in a nuclei? (1 to consider all nuclei and 2 to remove possible male cells)\n');

%% Extract data from file

ID = readtable([filenum '.csv']);

Allele = [];
NSNR = [];

% Reading the number of the last field from the long extension
LastIndex = strsplit(ID.File{length(ID.File)},'/');
Extension = strjoin(LastIndex(1:length(LastIndex)-1), '/');
LastIndex = strsplit(LastIndex{length(LastIndex)},'.');

for Fields = 1:str2double(LastIndex{1})
    
    if Fields < 10
        extension = [Extension, '/00', num2str(Fields), '.NM'];
    elseif Fields > 9 && Fields < 100
        extension = [Extension, '/0', num2str(Fields), '.NM'];
    else
        extension = [Extension, num2str(Fields), '.NM'];
    end
    
    PositionFile =  find(strcmp(extension, ID.File)); % Select current field of view

    nuclei = unique(ID.Nuclei(PositionFile)); % Identify the Nuclei to be analysed
    
    for n = 1:length(nuclei)
        PositionNuclei = PositionFile(eq(nuclei(n), ID.Nuclei(PositionFile))); % Select current nuclei
        
        preSNR = ID.NSNR(PositionNuclei(strcmp(ID.Channel(PositionNuclei), color(3))));
        [I, preSNR] = get_Intensity(ID.Channel(PositionNuclei), ID.Value(PositionNuclei), ...
            [ID.x(PositionNuclei)*0.13, ID.y(PositionNuclei)*0.13, ID.z(PositionNuclei)*0.3], Dye, color, MinDots, preSNR);
        
        NSNR = [NSNR; preSNR];
        
        Allele = [Allele; I];
    end
end

%% Calculate FRET efficiency

% Correct intensity contribution from donor and acceptor channels
if strcmp(Dye, 'B') && strcmpi(Correction, 'y')
    Allele(:, 3) = Allele(:, 3) - Allele(:, 2)*CorrDonor - Allele(:, 1)*CorrAcceptor;
end

% Calculate FRET ratio and efficiency
RatioFRETdonor = Allele(:, 3)./Allele(:, 2);
RatioFRETacceptor = Allele(:, 3)./Allele(:, 1);
EfficiencyFRET = Allele(:, 3)./(Allele(:, 3)+Allele(:, 2));


%% Plot figure

figure
subplot(2, 2, 1)
histogram(Allele(:, 1), 50)
xlabel('Acceptor Intensity')
legend(num2str(sum(~isnan(Allele(:, 1)))))
title(['median = ' num2str(round(nanmedian(Allele(:, 1)), 2))])

subplot(2, 2, 2)
histogram(Allele(:, 2), 50)
legend(num2str(sum(~isnan(Allele(:, 2)))))
xlabel('Donor Intensity')
title(['median = ' num2str(round(nanmedian(Allele(:, 2)), 2))])

subplot(2, 2, 3)
histogram(Allele(:, 3), 50)
legend(num2str(sum(~isnan(Allele(:, 3)))))
xlabel('FRET Intensity')
title(['median = ' num2str(round(nanmedian(Allele(:, 3)), 2))])

% Depends of which sample selected wether calculates efficiency or ratio
subplot(2, 2, 4)
if strcmp(Dye, 'D')
    histogram(RatioFRETdonor, 30)
    legend(num2str(sum(~isnan(RatioFRETdonor))))
    xlabel('FRET ratio with Donor')
    title(['median = ' num2str(round(nanmedian(RatioFRETdonor), 3))])
    
elseif strcmp(Dye, 'A')
    histogram(RatioFRETacceptor, 30)
    legend(num2str(sum(~isnan(RatioFRETacceptor))))
    xlabel('FRET ratio with Acceptor')
    title(['median = ' num2str(round(nanmedian(RatioFRETacceptor), 23))])

elseif strcmp(Dye, 'B')  
    histogram(EfficiencyFRET, 50)
    legend(num2str(sum(~isnan(EfficiencyFRET))))
    xlabel('FRET Efficiency')
    title(['median = ' num2str(round(nanmedian(EfficiencyFRET)*100, 1)) ' %'])
end

%% plot Ratio by no NSNR restriction or with <3 filter
if strcmp(Dye, 'A')
    figure
    subplot(1, 2, 1)
    scatter(NSNR, RatioFRETacceptor, 30, 'filled','MarkerFaceAlpha',0.6)
    xlabel('NSNR FRET')
    ylabel('FRET ratio with Acceptor')
    title(['median = ' num2str(nanmedian(RatioFRETacceptor))])
    grid
    
    subplot(1, 2, 2)
    scatter(NSNR(NSNR>3), RatioFRETacceptor(NSNR>3), 30, 'filled','MarkerFaceAlpha',0.6)
    xlabel('NSNR FRET')
    ylabel('FRET ratio with Acceptor')
    title(['median = ' num2str(nanmedian(RatioFRETacceptor(NSNR>3)))])
    legend(num2str(length(NSNR(NSNR>3))))
    grid
    
elseif strcmp(Dye, 'D')
    figure
    subplot(1, 2, 1)
    scatter(NSNR, RatioFRETdonor, 30, 'filled','MarkerFaceAlpha',0.6)
    xlabel('NSNR FRET')
    ylabel('FRET ratio with Acceptor')
    title(['median = ' num2str(nanmedian(RatioFRETdonor))])
    grid
    
    subplot(1, 2, 2)
    scatter(NSNR(NSNR>3), RatioFRETdonor(NSNR>3), 30, 'filled','MarkerFaceAlpha',0.6)
    xlabel('NSNR FRET')
    ylabel('FRET ratio with Donor')
    title(['median = ' num2str(nanmedian(RatioFRETdonor(NSNR>3)))])
    legend(num2str(length(NSNR(NSNR>3))))
    grid
end









function [Intensity, SNR] = get_Intensity(Channel, Values, position, Dye, color, MinDots, SNR)
% Remove dots positioned too close and consider only the brightest
% Remove values higher than 10^4
% Associate the correct dots from different channels
% Extract the intensities by channel
Intensity = nan(9, 3);

Values(Values > 10^4) =  nan;
Channel(Values > 10^4) =  {nan};
position(Values > 10^4, :) =  nan;

for i = 1:3
    
    if any(strcmp(Channel, color(i)))
        
        P = position(strcmp(Channel, color(i)), :);
        Distance = triu(squareform(pdist(P)));    
        [a, b] = find(Distance < 0.0 & Distance ~= 0); % 2µm as the threshold
        
        if ~isempty(a)
            %find minimum Intensity value and relate with current channel position then to the general positions
            vec = [a, b];
            Int = Values(strcmp(Channel, color(i)));
            toRemove = ismember(position, P(ismember(Int, min(Int(vec'))), :));
            position(toRemove) = nan; 
            Channel(toRemove(:, 1)) = {nan}; 
            Values(toRemove(:, 1)) = nan; 
        end
        
        Intensity(1:sum(strcmp(Channel, color(i))), i) = Values(strcmp(Channel, color(i)));
    end
end

% superior to one means that excludes nuclei with only 1 dot (male cells)
if sum(~isnan(Intensity(:, 3))) >= MinDots % It's not enough, should be a way to unconsider one unpaired case
    
    if strcmp(Dye, 'A') && sum(~isnan(Intensity(:, 1))) >= MinDots && (sum(strcmp(Channel, color(1))) == sum(strcmp(Channel, color(3)))) 
        Intensity(:, 2) = nan;
        
        Intensity(all(isnan(Intensity), 2), :) = [];
        
        [~, i] = sort(sum(position(strcmp(Channel, color(1)), :), 2));
        Intensity(:, 1) = Intensity(i, 1);
        [~, i] = sort(sum(position(strcmp(Channel, color(3)), :), 2));
        Intensity(:, 3) = Intensity(i, 3);
    
        SNR = SNR(i);
        
    elseif strcmp(Dye, 'D') && sum(~isnan(Intensity(:, 2))) >= MinDots  && (sum(strcmp(Channel, color(2))) == sum(strcmp(Channel, color(3))))
        Intensity(:, 1) = nan;
        
        Intensity(all(isnan(Intensity), 2), :) = [];
        
        [~, i] = sort(sum(position(strcmp(Channel, color(2)), :), 2));
        Intensity(:, 2) = Intensity(i, 2);
        [~, i] = sort(sum(position(strcmp(Channel, color(3)), :), 2));
        Intensity(:, 3) = Intensity(i, 3);
        
        SNR = SNR(i);
        
    elseif strcmp(Dye, 'B') && sum(~isnan(Intensity(:, 2))) >= MinDots && sum(~isnan(Intensity(:, 1))) >= MinDots  && (sum(strcmp(Channel, color(2))) == sum(strcmp(Channel, color(3)))) && (sum(strcmp(Channel, color(1))) == sum(strcmp(Channel, color(3))))
        
        Intensity(all(isnan(Intensity), 2), :) = [];
        
        [~, i] = sort(sum(position(strcmp(Channel, color(1)), :), 2));
        Intensity(:, 1) = Intensity(i, 1);
        [~, i] = sort(sum(position(strcmp(Channel, color(2)), :), 2));
        Intensity(:, 2) = Intensity(i, 2);
        [~, i] = sort(sum(position(strcmp(Channel, color(3)), :), 2));
        Intensity(:, 3) = Intensity(i, 3);
        
        SNR = SNR(i);
    
    else
        Intensity = nan(9, 3);
        SNR = nan(9, 1);
    end
    
else
    Intensity = nan(9, 3);
    SNR = nan(9, 1);
end


if any(any(Intensity==0))
    disp('problem')
end
end
