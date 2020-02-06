

% APTで撮影したファイルのフォルダーを指定、フォルダー直下にB_*、D_*、L_* ファイル
APTFolder = '/Users/tjm/Pictures/PixInsight/APT_Images/2019-12-31/MATLABTEST/';
% APTFolder = uigetdir();
% APTのfitsはdouble

% BIAS = reading noise by Gain and Offset
% DARK = Light Length Thermal noise + BIAS
% LIGHT = (FlatDistortion Signal) + Light Length Thermal noise + BIAS

% FLAT = FlatDistortion + FLAT Length Thermal noise + BIAS
% FlatDistortion = FLAT - BIAS - Light Length Thermal noise / ( Light Length / FLAT Length )

% フォルダーの作成
% MASTER BIAS DARK FLAT LIGHT のマスターファイルを保存します
[~, ~, ~] = mkdir(APTFolder, 'MASTER');
dateformat = 'yyyy/mm/dd HH:MM:SS.FFF';
fprintf('%s %s\n',datestr(now,dateformat), 'process started'); 

dlist = dir(APTFolder);
% remove . and .. 
dlist = dlist(cellfun(@(x) x ~= 1, {dlist.isdir}));

BiasStack = [];
for i = 1:length(dlist)
    if contains(dlist(i).name, 'B_')
        fprintf('%s %s %s\n',datestr(now,dateformat), 'reading BIAS', dlist(i).name);
        dt = fitsread(fullfile(APTFolder, dlist(i).name));
        BiasStack = cat(3, BiasStack, dt);
    end
end

fprintf('%s %s\n', datestr(now,dateformat),'Master BIAS Start');
% (2) Stack Average, Sigma Clipping (-3σ, +3σ)
B_Ave = mean(BiasStack, 3);
B_Mid = median(BiasStack, 3);
B_Std = std(BiasStack, 0, 3);

% Rejection Map, 1 = Reject / 0 = Keep
Rejection = (BiasStack < B_Ave - B_Std .* 3 |...
             B_Ave + B_Std .* 3 < BiasStack);

% Keep Original + Reject and Replace by Median
BiasStack = BiasStack .* ~Rejection + B_Mid .* Rejection;    

BiasOutput = mean(BiasStack, 3);
fitswrite(BiasOutput, fullfile(APTFolder, 'MASTER', 'Master_Bias.fit'));
fprintf('%s %s\n', datestr(now,dateformat),'Master BIAS done');

% Dark FIle, same as BIAS
DarkStack = [];
for i = 1:length(dlist)
    if contains(dlist(i).name, 'D_')
        fprintf('%s %s %s\n',datestr(now,dateformat), 'reading DARK', dlist(i).name); 
        dt = fitsread(fullfile(APTFolder, dlist(i).name));
        DarkStack = cat(3, DarkStack, dt);
    end
end

fprintf('%s %s\n', datestr(now,dateformat),'Master DARK Start');
% (2) Stack Average, Sigma Clipping (-3σ, +3σ)
D_Ave = mean(DarkStack, 3);
D_Mid = median(DarkStack, 3);
D_Std = std(DarkStack, 0, 3);

% Rejection Map, 1 = Reject / 0 = Keep
Rejection = (DarkStack < D_Ave - D_Std .* 3 |...
             D_Ave + D_Std .* 3 < DarkStack);

% Keep Original + Reject and Replace by Median
DarkStack = DarkStack .* ~Rejection + D_Mid .* Rejection;    

DarkOutput = mean(DarkStack, 3);
fitswrite(DarkOutput, fullfile(APTFolder, 'MASTER', 'Master_Dark.fit'));
fprintf('%s %s\n', datestr(now,dateformat),'Master DARK done');

LightStack = [];
for i = 1:length(dlist)
    if contains(dlist(i).name, 'L_')
        fprintf('%s %s %s\n',datestr(now,dateformat), 'reading LIGHT', dlist(i).name); 
        
        dt = fitsread(fullfile(APTFolder, dlist(i).name));
        % Light File calibration, Light = Light - Dark

        % DarkOutputはオフセット1976のDouble値
        % 読み込んだdtはだいたいオフセット1977の整数値でフォーマットはDouble
        dtAvg = mean(mean(dt));
        dtStd = std(std(dt));
        dt = dt - DarkOutput;
        % 負の値を置き換える
        Negatives = (dt < (dtAvg - 5 * dtStd));
        dt = dt .* ~Negatives + dtAvg .* Negatives;
        % 0-65535に正規化
        dtMin = min(min(dt));
        dtMax = max(max(dt));
        dt = (dt - dtMin) .* 65535 ./ (dtMax - dtMin) ;

        % Light File Debayer, rggb (int16)
        dt = cast(dt, 'uint16');
        im2 = demosaic(dt, 'rggb');  
        % Now Light Flame is NxNx3
        LightStack = cat(4, LightStack, im2);
    end
end

fprintf('%s %s\n', datestr(now,dateformat),'Light Alignment');
[~,~,~,len] = size(LightStack);
% Light File alignment
imFixed = LightStack(:,:,:,round(len / 2)); % reference image
imFixed = imadjust(imFixed, stretchlim(imFixed, [0.0 1.0]));
Areafixed = imref2d(size(imFixed)); % reference coordinates
AlignedStack = [];
for i = 1:len
    imMoving = LightStack(:,:,:,i);
    imMoving = imadjust(imMoving, stretchlim(imMoving, [0.0 1.0]));
    tFormation = imregcorr(imMoving,imFixed,'translation');
    fprintf('%s %s %i\n', datestr(now,dateformat),'Light transform Frame:', i);
    disp(tFormation.T);
    AlignedStack = cat(4, AlignedStack, imwarp(imMoving,tFormation,'OutputView',Areafixed));
end    
AlignedStack = cast(AlignedStack, 'single');
fprintf('%s %s\n', datestr(now,dateformat),'Light Stacking');
% (2) Stack Average, Sigma Clipping (-4σ, +3σ)
L_Ave = mean(AlignedStack, 4);
L_Mid = median(AlignedStack, 4);
L_Std = std(AlignedStack, 0, 4);
% Rejection Map, 1 = Reject / 0 = Keep
Rejection = (AlignedStack < L_Ave - L_Std .* 4 |...
             L_Ave + L_Std .* 3 < AlignedStack);
% Keep Original + Reject and Replace by Median
AlignedStack = AlignedStack .* ~Rejection + L_Mid .* Rejection;  
% Light File Output
LightOutput = mean(AlignedStack, 4);

fitswrite(LightOutput, fullfile(APTFolder, 'MASTER', 'Master_LIGHT.fit'));
fprintf('%s %s\n', datestr(now,dateformat),'Master LIGHT done');

fprintf('%s %s\n', datestr(now,dateformat),'process end');
