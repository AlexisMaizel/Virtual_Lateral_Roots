function cellStartIndex = determineDeformationDirection( startEnd, lineX, lineY )
realS = startEnd(1, 1:2);
realE = startEnd(1, 4:5);
defS = [ lineX(1) lineY(1) ];
defE = [ lineX(2) lineY(2) ];

% nuclei deformation
realDir = realE - realS;
realDir = normalize(realDir);
% possible deformation directions
defDir1 = defE - defS;
defDir2 = defS - defE;
defDir1 = normalize(defDir1);
defDir2 = normalize(defDir2);

angle1 = acos( dot(realDir, defDir1) );
angle2 = acos( dot(realDir, defDir2) );

if angle1 < angle2
  cellStartIndex = 1;
else
  cellStartIndex = 2;
end

