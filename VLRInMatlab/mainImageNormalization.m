setWorkingPathProperties()

figure;
plotIndex = 1;

for t=3:23
  if t < 10
    digit = '00';
  elseif t < 100
    digit = '0';
  else
    digit = '';
  end
  nucleiFileName = strcat( 'I:\NewDatasets\Zeiss\20160426\red\cropped_spim_TL', digit, num2str(t), '_Angle1.tif' );
  %membraneFileName = strcat( 'I:\NewDatasets\Zeiss\20160426\green\cropped_spim_TL006_Angle1.tif' );
  
  % create 3D array of binary image data
  imageStack = readTIFstack( char(nucleiFileName) );
  height = size( imageStack, 1 );
  width = size( imageStack, 2 );
  slices = size( imageStack, 3 );
  
  if t == 6 || t == 16 || t == 17 || t == 18 || t == 19 || t == 23
    nImageStack = imageStack;
    
    % determine min and max values
    minI = min( nImageStack(:) );
    maxI = max( nImageStack(:) );
    percVal = 0.2;
    barrierI = maxI * ( 1. - percVal );
    percAdd = 0.1;
    addI = percAdd * maxI;
    
    for s=1:slices
      for j=1:width
        for i=1:height
          val = nImageStack( i, j, s );
          if val > barrierI
            nImageStack( i, j, s ) = val - addI;
          else
            nImageStack( i, j, s ) = val + addI;
          end
        end
      end
    end
    
    imageDir = strcat( 'I:\NewDatasets\Zeiss\20160426\red\partiallyNormalized\cropped_spim_TL', digit, num2str(t), '.tif' );
    %nImageStack = uint16(nImageStack).*65535;
    options.message = true;
    writeTIFstack( nImageStack, char(imageDir), 2^31, options );
  else
    imageDir = strcat( 'I:\NewDatasets\Zeiss\20160426\red\partiallyNormalized\cropped_spim_TL', digit, num2str(t), '.tif' );
    options.message = true;
    writeTIFstack( imageStack, char(imageDir), 2^31, options );
  end
  
%   subplot( 2, 2, plotIndex )
%   h = showMIP( imageStack );
%   plotIndex = plotIndex + 1;
%   subplot( 2, 2, plotIndex )
%   h2 = showMIP( nImageStack );
%   plotIndex = plotIndex + 1; 
end