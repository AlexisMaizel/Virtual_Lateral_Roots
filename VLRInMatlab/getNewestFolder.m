function newestFolder = getNewestFolder( folderPath )
allFolders = dir( folderPath );
dirFlags = [ allFolders.isdir ];
allFolders = allFolders( dirFlags );
newestFolder = {};
timeDiff = 1000;
for f=1:size( allFolders, 1 )
  name = allFolders( f, 1 ).name;
  fullDate = allFolders( f, 1).date;
  diff = datenum( datetime('now') ) - datenum( fullDate );
  if diff < timeDiff
    timeDiff = diff;
    newestFolder = name;
  end
end