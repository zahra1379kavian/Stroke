%%Phase Locking Value(PLV)
clc; clear; %close all;

data = cell2mat(struct2cell(load('sub02')));
data = double(data);

filtSpec.order = 50;
filtSpec.range = [8 12]; %alpha range
plv = PLV_func(data, filtSpec);
save PLV_alpha plv

filtSpec.range = [12 30]; %beta range
plv = PLV_func(data, filtSpec);
save PLV_beta plv

filtSpec.range = [30 100]; %gamma range
plv = PLV_func(data, filtSpec);
save PLV_gamma plv
%%
data = cell2mat(struct2cell(load('PLVSub21')));
meanplv = squeeze(mean(data,1));

figure
imagesc(meanplv);
xlabel('Electrodes')
ylabel('Electrodes')
title('Subject 26')
colormap(parula(3))


%%
%%%XY coordinate
elocx = cell2mat(struct2cell(load('X_channel_location')));
elocy = cell2mat(struct2cell(load('Y_channel_location')));

plotplv_scatter(meanplv,0.7,elocx,elocy,[1 0 0])
plotplv_scatter(meanplv,0.8,elocx,elocy,[0 1 0])
plotplv_scatter(meanplv,0.9,elocx,elocy,[0 0 1])


function [plv] = PLV_func(eegData, filtSpec)
numChannels = size(eegData, 1);
trial = size(eegData, 2);
numTrials = size(eegData, 3);

bpFilt = designfilt('bandpassfir','FilterOrder',filtSpec.order, ...
         'CutoffFrequency1',filtSpec.range(1),'CutoffFrequency2',filtSpec.range(2), ...
         'SampleRate',1000);
filteredData = zeros(size(eegData));
for i = 1:numChannels
filteredData(i,:,:) = filtfilt(bpFilt,squeeze(eegData(i,:,:)));
end
%filteredData = eegData;

angle_filteredData = zeros(size(eegData));
for channelCount = 1:numChannels
    angle_filteredData(channelCount, :, :) = angle(hilbert(filteredData(channelCount, :, :)));
end

plv = zeros(numTrials, numChannels, numChannels);

for channelCount = 1:numChannels
    for compareChannelCount = 1:numChannels
            if channelCount ~= compareChannelCount
            plv(:, channelCount, compareChannelCount) = abs(mean(squeeze(exp(1i*(angle_filteredData(channelCount,:,:) - angle_filteredData(compareChannelCount,:,:)))), 1));
            end
    end
end
end


function plotplv_scatter(meanplv,threshold1,elocx,elocy,C)
s = 1;
s1 = 1;

for i = 1:size(meanplv,1)
    for j = 1:size(meanplv,2)
        if meanplv(i,j)>=threshold1 & meanplv(i,j)~=1
            x(s) = i;
            y(s) = j;
            v(s) = meanplv(i,j);
            s = s+1;
        end
    end
end

figure
scatter(elocx,elocy,100,[206/255 206/255 206/255],'filled')

hold on
for i = 1:size(x,2)   
    line([elocx(x(1,i),1),elocx(y(1,i),1)],[elocy(x(1,i),1),elocy(y(1,i),1)],'linewidth',v(1,i)+2,'Color',C);
    hold on
end

end

