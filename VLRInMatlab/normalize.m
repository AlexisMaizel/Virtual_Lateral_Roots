function v = normalize(vec)
dimR = size(vec, 1);
dimC = size(vec, 2);
if dimR > dimC
  dim = dimR;
else
  dim = dimC;
end
norm = 0;
for i=1:dim
  if dim == dimR
    norm = norm + vec(i,:).*vec(i,:);
  else
    norm = norm + vec(:,i).*vec(:,i);
  end
end

norm = sqrt(norm);

v = vec./norm;