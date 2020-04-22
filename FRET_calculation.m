clear all
close all
clc

%% Input information
% Check the correction factors and field resolution

filenum = input('File name? (e.g. iAM382.csv)\n', 's');
Dye = input('Type of dyes? (e.g. Acceptor or Donor or Both)\n', 's');

color = input('Which dyes? (C for Cy3-Cy5 and A for a488-a594)\n', 's');

if strcmpi(color, 'C')
    color={'Cy5'; 'tmr'; 'tmrfret'};
elseif strcmpi(color, 'A')
    color={'a594'; 'a488'; 'a594fret'};
end

Correction = input('Correct images? (y/n)\n', 's');

if strcmpi(Correction, 'y')
    CorrDonor = input('Correction ratio for Donor intensity? (e.g. for Cy3 is 0.278 and a488 is 0.085)\n');
    CorrAcceptor = input('Correction ratio for Acceptor intensity? (e.g. for Cy5 is 0.09 and a594 is 0.056)\n');
end

disp('If there is a group of dots at a distance lower than 2 µm, only the brigher dot is considered.')
disp('Dots with aberrant intensity values are removed. (e.g. 10000 when normal is 10)')
MinDots = input('Minimum dots in a nuclei? (1 to consider all nuclei and 2 to remove possible male cells)\n');

%% Extract data from file

ID = readtable(filenum);

Allele = [];

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
        
        Allele = [Allele; get_Intensity(ID.Channel(PositionNuclei), ID.Value(PositionNuclei), ...
            [ID.x(PositionNuclei)*0.13, ID.y(PositionNuclei)*0.13, ID.z(PositionNuclei)*0.3], Dye, color, MinDots)];
    end
end

%% Calculate FRET efficiency

% Correct intensity contribution from donor and acceptor channels
if strcmp(Dye, 'Both') && strcmpi(Correction, 'y')
    Allele(:, 3) = Allele(:, 3) - Allele(:, 2)*CorrDonor - Allele(:, 1)*CorrAcceptor;
end

% Calculate FRET ratio and efficiency
RatioFRETdonor = Allele(:, 3)./Allele(:, 2);
RatioFRETacceptor = Allele(:, 3)./Allele(:, 1);
EfficiencyFRET = Allele(:, 3)./(Allele(:, 3)+Allele(:, 2));


%% Plot figure

figure
subplot(2, 2, 1)
histogram(Allele(:, 1), 20)
xlabel('Acceptor Intensity')
legend(num2str(sum(~isnan(Allele(:, 1)))))
title(['median = ' num2str(round(nanmedian(Allele(:, 1)), 2))])

subplot(2, 2, 2)
histogram(Allele(:, 2), 20)
legend(num2str(sum(~isnan(Allele(:, 2)))))
xlabel('Donor Intensity')
title(['median = ' num2str(round(nanmedian(Allele(:, 2)), 2))])

subplot(2, 2, 3)
histogram(Allele(:, 3), 20)
legend(num2str(sum(~isnan(Allele(:, 3)))))
xlabel('FRET Intensity')
title(['median = ' num2str(round(nanmedian(Allele(:, 3)), 2))])

subplot(2, 2, 4)
% Only consider the ratio comparison with the channel active
if strcmp(Dye, 'Donor')
    histogram(RatioFRETdonor, 30)
    legend(num2str(sum(~isnan(RatioFRETdonor))))
    xlabel('FRET ratio with Donor')
    title(['median = ' num2str(round(nanmedian(RatioFRETacceptor), 2))])
    
elseif strcmp(Dye, 'Acceptor')
    histogram(RatioFRETacceptor, 30)
    legend(num2str(sum(~isnan(RatioFRETacceptor))))
    xlabel('FRET ratio with Acceptor')
    title(['median = ' num2str(round(nanmedian(RatioFRETacceptor), 2))])

elseif strcmp(Dye, 'Both')  
    histogram(EfficiencyFRET, 50)
    legend(num2str(sum(~isnan(EfficiencyFRET))))
    xlabel('FRET Efficiency')
    title(['median = ' num2str(round(nanmedian(EfficiencyFRET)*100, 1)) ' %'])
end













function Intensity = get_Intensity(Channel, Values, position, Dye, color, MinDots)
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
        [a, b] = find(Distance < 2.0 & Distance ~= 0); % 2µm as the threshold
        
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
    
    if strcmp(Dye, 'Acceptor') && sum(~isnan(Intensity(:, 1))) >= MinDots  
        Intensity(:, 2) = nan;
        
        Intensity = closestDot(position(strcmp(Channel, color(1)), :), ...
            position(strcmp(Channel, color(3)), :), Intensity, 1);
    
        
    elseif strcmp(Dye, 'Donor') && sum(~isnan(Intensity(:, 2))) >= MinDots
        Intensity(:, 1) = nan;
        
        Intensity = closestDot(position(strcmp(Channel, color(2)), :), ...
            position(strcmp(Channel, color(3)), :), Intensity, 2);
        
        
    elseif strcmp(Dye, 'Both') && sum(~isnan(Intensity(:, 2))) >= MinDots && sum(~isnan(Intensity(:, 1))) >= MinDots
        
        Intensity = closestDot(position(strcmp(Channel, color(2)), :), ...
            position(strcmp(Channel, color(3)), :), Intensity, 2);
        
        Intensity = closestDot(position(strcmp(Channel, color(1)), :), ...
            position(strcmp(Channel, color(3)), :), Intensity, 1);
    
    else
        Intensity = nan(9, 3);
    end
    
else
    Intensity = nan(9, 3);
end


if any(any(Intensity==0))
    disp('problem')
end
end



 function Intensity = closestDot(AorD, FRET, Intensity, num)
        AorD = sum(AorD, 2);
        FRET = sum(FRET, 2);
        
        A = repmat(AorD, 1, length(FRET));
        [~,closestIndex] = min(abs(A-FRET'));
        
        Intensity(closestIndex, num)=Intensity(1:length(closestIndex), num);
 end