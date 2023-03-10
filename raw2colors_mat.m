function outmat = raw2colors_mat(regularRaw,colorChange)

% regularRaw = regularized rawassn matrix from link Infomap
% colorChange = matrix of color changes (formats is col 1 change to col 2)
regularNodes = regularRaw;

dataColor = regularNodes;
for i=1:size(regularNodes,2)
    temp = dataColor(:,i);
    for j=1:size(colorChange,1)
        inds = find(regularNodes(:,i)==colorChange(j,1));
        temp(inds)=colorChange(j,2);
    end
    dataColor(:,i)=temp; 
end

outmat = dataColor;
