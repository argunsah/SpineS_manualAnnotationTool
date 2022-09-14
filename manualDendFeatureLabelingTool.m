close all
clear all
clc

currentDir = pwd;
addpath(genpath(currentDir));

desiredPixelSize    = 0.035;        % Do not change this value! 
fullFileList        = uipickfiles;

labName             = 'turtleLab';  % Change this to your Lab Name

pixelSizeZ          = 0.3;          % Change this for tif files because tif doesn't keep this value in it. 
                                    % For nd2 no change necessary, it will
                                    % be overwritten by line 28

for fileNo = 1:size(fullFileList,2)
    
    [dataInfo, multiCube]       = readMicroscopyData(fullFileList(fileNo));
    [filePath,fileName,fileExt] = fileparts(fullFileList{fileNo});
    saveFileName                = sprintf('%s_%s',fileName,'Annotations');
    
    try
        pixelSizeXY = 1/(dataInfo.XResolution);
    catch
        try
            pixelSizeXY = (dataInfo.voxSizeX)*10^6;
        catch
            try
                pixelSizeXY = dataInfo.ScaleX;
                pixelSizeZ  = dataInfo.ScaleZ;
            catch
                pixelSizeXY = nan;
            end
        end
    end
    
    % Data Specs
    convFactor          = 1/(desiredPixelSize/pixelSizeXY);
    
    % Extract the Cube XxYxZ
    timePoint           = 1;
    imageChannel        = 1;
    imageCube           = squeeze(multiCube(:,:,:,imageChannel,timePoint)); %XYZCT
    imageCube           = imageCube-min(imageCube(:));
    temp                = double(imageCube);
    temp                = temp/max(temp(:));
    
    imageCube16btScaled = uint16(temp*2^16-1);
    
    [mipIm,indIm]       = max(imageCube16btScaled,[],3);
    tempMIP             = imresize(mipIm,convFactor);
    tempZ               = imresize(indIm,convFactor);
    tempZnorm           = tempZ/size(imageCube,3);

    figure, imshow(imadjust(tempMIP),[]);
    h = imcontrast(gca);
    waitfor(h);

    inputOptions = {'Dendrite Tracing','Annotation'};
    defSelection = inputOptions{2};
    iSel = bttnChoiseDialog(inputOptions, '', defSelection,'What do you want to do?');

    if iSel==1
        % Tracing Section
        clear x;
        clear y;
    
        b      = 10000;
        xx     = 1;
        while xx < b
            clear x;
            clear y;
            clear z;
           
            [x,y] = ginputzoom(size(tempMIP,2),size(tempMIP,1));
            
            if length(x)>1
           
                hA = gca;
                resetplotview(hA,'InitializeCurrentView');
                set(hA,'xlim',[1 size(tempMIP,2)]);
                set(hA,'ylim',[1 size(tempMIP,1)]);
    
                for k = 1:length(x)
                    z(k) = tempZ(round(y(k)),round(x(k)));
                end
    
                meanSL = [];
                for n = 1:length(x)-1
                    SourcePoint    = [round(y(n))   round(x(n))   z(n)]';
                    StartPoint     = [round(y(n+1)) round(x(n+1)) z(n+1)]';
    
                    T1_FMM1        = msfm2d((1000000*double(tempMIP)+1000), SourcePoint(1:2), true, true);
                    T1_FMM2        = msfm2d((1000000*double(tempMIP)+1000), StartPoint(1:2) , true, true);
    
                    ShortestLine   = shortestpath(T1_FMM1,StartPoint(1:2),SourcePoint(1:2));
                    ShortestLine2  = shortestpath(T1_FMM2,SourcePoint(1:2),StartPoint(1:2));
    
                    meanSL         = [meanSL ;meanOfTwoLines(ShortestLine,ShortestLine2)];
                end

                subTemp = sub2ind(size(tempMIP),meanSL(:,1),meanSL(:,2));
                allZ    = improfile(tempZ,meanSL(:,2),meanSL(:,1),size(meanSL,1));
                allZ    = sgolayfilt(allZ,3,2*round(size(meanSL,1)/6)-1);
                hh      = scatter(meanSL(:,2),meanSL(:,1),1,'.r');
                drawnow;

                allZ_microMeter     = allZ*pixelSizeZ;
                meanSL_microMeter   = meanSL*desiredPixelSize;

                traceProfile = [];

                traceProfile(:,1) = meanSL_microMeter(:,2);
                traceProfile(:,2) = meanSL_microMeter(:,1);
                traceProfile(:,3) = allZ_microMeter;

                D   = diff(traceProfile,1,1);
                processLength_microMeter = double(trace(sqrt(D*D.')));
    
                answer = questdlg('Do you like the Tracing?', ...
                    'Tracing Quality', ...
                    'Yes','No','Done','Yes');
                % Handle response
                switch answer
                    case 'Yes'
                        branchLen{1,xx}         = processLength_microMeter;
                        numDendrites            = xx;
                        shortestLinesAll{1,xx}  = traceProfile;
                        xx = xx + 1;
                    case 'No'
                        delete(hh);
                    case 'Done'
                        branchLen{1,xx}         = processLength_microMeter;
                        numDendrites            = xx;
                        shortestLinesAll{1,xx}  = traceProfile;
   
                        b = xx;

                        savefile = sprintf('Tracing-%s-%s-rnd%d.mat',labName,saveFileName,round(rand(1)*1000));
                        save(savefile,'tempMIP','tempZ','convFactor',...
                            'branchLen','numDendrites','pixelSizeZ',...
                            'shortestLinesAll','desiredPixelSize','pixelSizeXY');
                end
            end
        end
    else
        % Annotation Section
        iSel2   = 0;
        x_peaks = [];
        y_peaks = [];
        trainingLabels = [];

        while iSel2<10
            labelList = {'Single Spine','Multi Spine',...
                'Neck Base','Dendrite Edge','Dendrite',...
                'Spine Neck','Spine Head Edge','Bouton or Axon','Noise','Done!'};
            defSelection2 = labelList{10};

            iSel2 = bttnChoiseDialog(labelList,...
                '', defSelection2,'What do you want to annotate?');

            if iSel2==10
                padValue            = 49;
                mipIm_padded        = padarray(tempMIP ,[padValue padValue],0,'both');
                zIm_padded          = padarray(tempZ   ,[padValue padValue],0,'both');
            
                xp_padded           = x_peaks + padValue;
                yp_padded           = y_peaks + padValue;
                
                cubes               = cell2mat(...
                    arrayfun(@(x,y) reshape(mipIm_padded((x-(padValue-1)):(x+(padValue-1)),...
                    (y-(padValue-1)):(y+(padValue-1))),...
                    [padValue*2-1,padValue*2-1]),...
                    round(xp_padded),round(yp_padded),'uniformoutput',false));
                cubes               = permute(cubes,[2 1]);
                cubes               = reshape(cubes', [padValue*2-1 padValue*2-1 length(x_peaks)]);
            
                cubesZ              = cell2mat(...
                    arrayfun(@(x,y) reshape(zIm_padded((x-(padValue-1)):(x+(padValue-1)),...
                    (y-(padValue-1)):(y+(padValue-1))),...
                    [padValue*2-1,padValue*2-1]),...
                    round(xp_padded),round(yp_padded),'uniformoutput',false));
                cubesZ              = permute(cubesZ,[2 1]);
                cubesZ              = reshape(cubesZ', [padValue*2-1 padValue*2-1 length(x_peaks)]);
            
                savefile = sprintf('TrainingSet-%s-%s-rnd%d.mat',labName,saveFileName,round(rand(1)*1000));
                save(savefile,'tempMIP','tempZ','convFactor',...
                    'padValue','x_peaks','y_peaks',...
                    'xp_padded','yp_padded','cubes',...
                    'cubesZ','trainingLabels','labelList',...
                    'desiredPixelSize','pixelSizeXY','pixelSizeZ');
            else
                clear X;
                clear Y;
            
                title(sprintf('Annotate %s',labelList{iSel2}));
                [Y,X]   = ginputzoom(size(tempMIP,2),size(tempMIP,1)); 
                x_peaks = [x_peaks round(X)];
                y_peaks = [y_peaks round(Y)];
               
                hA = gca;
                resetplotview(hA,'InitializeCurrentView');
                set(hA,'xlim',[1 size(tempMIP,2)]);
                set(hA,'ylim',[1 size(tempMIP,1)]);
    
                trainingLabels = [trainingLabels iSel2*ones(1,length(X))];
            end
        end
    end
end