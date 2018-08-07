close all;
clear all;

rawFolder = '/home/mobvoi/work/data/hixiaowen_lantuo_spk_1/';

rawPath = sprintf('%s/*.wav',rawFolder);
file = dir(rawPath);
fileNum = size(file, 1);

lpowers = zeros(fileNum, 1);

for ii=1:fileNum
    filename = file(ii).name;
    filename = strcat(rawFolder, filename);
    str = sprintf('num %d %s process started\n',ii, filename);
    [a, b] = vad(filename);
    lpowers(ii) = b;
end;

%%%%%%%%%%%%%%%

rawFolder = '/home/mobvoi/work/data/hixiaowen_mobvoi_spk_1/';

rawPath = sprintf('%s/*.wav',rawFolder);
file = dir(rawPath);
fileNum = size(file, 1);

mpowers = zeros(fileNum, 1);

for ii=1:fileNum
    filename = file(ii).name;
    filename = strcat(rawFolder, filename);
    str = sprintf('num %d %s process started\n',ii, filename);
    [a, b] = vad(filename);
    mpowers(ii) = b;
end;

figure(1);
plot(lpowers, 'r'); hold on;
plot(mpowers, 'g');