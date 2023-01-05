function [dataInfo, multiCube] = readMicroscopyData(fullFileList)

    % Check if it is a list of Files or Folders
    % Data Order : XYCZT
    numTimePoints = size(fullFileList,2);

    for numTime = 1:numTimePoints

        dataName     = fullFileList{numTime};
        [~, ~, fExt] = fileparts(dataName);

        if ~isempty(fExt)
            switch lower(fExt)
                case '.lsm'
                    dataInfo = lsmread(dataName,'InfoOnly');
                    data     = lsmread(dataName); %[T(ime), C(hannel), Z, X, Y]
                    data     = permute(data,[4 5 3 2 1]);
                case {'.czi','.oib','.oif','.nd2'}
                    out      = ReadImage6D(dataName, true, 1); %XYCZT
                    dataInfo = out{2};
                    data     = out{1};
                case {'.tif','.tiff'}
                    % tagNames = tiffObj.getTagNames()
                    tiffObj  = Tiff(dataName,'r');
                    data     = FastTiff(dataName);
                    try
                        dataInfo.XResolution = tiffObj.getTag('XResolution');
                    catch
                        dataInfo.XResolution = nan;
                    end
            end
        else
            filelist            = dir(fullfile(dataName,'*.tif*'));
            
            reader              = bfGetReader(fullfile(dataName,filelist(1).name));
            omeMetadata         = reader.getMetadataStore;
            voxelSizeX          = omeMetadata.getPixelsPhysicalSizeX(0).value; % in ?m
            voxelSizeXdouble    = voxelSizeX.doubleValue(); 
            dataInfo.ScaleX     = voxelSizeXdouble;

            for z = 1:size(filelist,1)
                temp = FastTiff(fullfile(dataName,filelist(z).name));
                data(:,:,z) = temp;
            end
        end

        if numTime == 1 && numTimePoints>1
            multiCube = zeros(size(data,1),size(data,2),dataInfo.SizeC,dataInfo.SizeZ,numTimePoints); %XYZCT
        end
        multiCube(:,:,:,:,numTime) = squeeze(data);
    end
end
