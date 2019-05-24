clear all
close all
clc

NPC129 = readtable('GSM1828646_NPC_wt_129_merge.txt', 'Format','%s%f%f%f'); % Inactive allele format is chromosome, start, end, reads normalized
NPCCast = readtable('GSM1828646_NPC_wt_Cast_merge.txt', 'Format','%s%f%f%f'); % Active allele

ChrX129 = NPC129(contains(NPC129{:, 1}, 'chrX'), :); % Extract the information corresponding to chrX only
Xist129 = find((ChrX129{:, 2} > 98*10^6) & (ChrX129{:, 2} < 105*10^6)); % Get the position of the information in between 98 and 105 Mbp
Xi = ChrX129{Xist129, 4}; % Only the reads

ChrXCast = NPCCast(contains(NPCCast{:, 1}, 'chrX'), :);
XistCast = find((ChrXCast{:, 2} > 98*10^6) & (ChrXCast{:, 2} < 105*10^6));
Xa = ChrXCast{XistCast, 4};

% 129 is the inactive allele so the reads longer since the DNA is less
% digested because has lack of exposure to enzymes due to the repression.
% Therefore, for a compehensive comparison the reads of Cast are included
% in the reads of 129.
Start129 = ChrX129{Xist129, 2}; % Only the 129 allele start position
StartCast = ChrXCast{XistCast, 2}; % Only the Cast allele start position
A = repmat(StartCast,[1 length(Start129)]); % Repeat in a matrix all the values from each allele
[minValue,closestIndex] = min(abs(A-Start129')); % find the closest value in between 129 and Cast
closestValue = StartCast(closestIndex);

% Only the first read of Xa in the region of Xi is considered. For all the
% Xa reads to be considered, the Xi value has to be repeated has many times
% the Xa was digested.
ratio = Xi./Xa(closestIndex);

Inter = ChrX129{Xist129, 2}-ChrX129{Xist129-1, 2}; % Space in between the start positions
F = NaN(length(ratio), 600); % Pre-allocation matrix
for Start_position = 1:length(ratio)
    Accumulative_distance = 0;
    Extension_distance = 1;
    Ratio_temp = [];
    
    while Accumulative_distance < 10^5 % Window value
        
        if Start_position + Extension_distance > length(ratio)
            break
        end
        
        Accumulative_distance = Accumulative_distance + Inter(Start_position+Extension_distance);
        Ratio_temp(Extension_distance) = ratio(Start_position+Extension_distance);
        Extension_distance = Extension_distance+1;
    end
    
    %Instead of the sum to give 500, it would be nicer to select a
    %threshold of actual peaks like less than 0.3 and higher than 3
    if Extension_distance ~= 1 && sum(Ratio_temp>3) > 16 && Accumulative_distance < 1.1*10^5 && sum(Ratio_temp<0.40) == 0
        F(Start_position, 1:length(Ratio_temp))=Ratio_temp;
        disp('window length')
        disp(Accumulative_distance);
        disp('mean ratio')
        disp(nanmean(Ratio_temp));
        disp('start window coordinate')
        disp(ChrX129{Xist129(Start_position), 2})
        figure(Start_position)
        hist(Ratio_temp, 30)
    end
end