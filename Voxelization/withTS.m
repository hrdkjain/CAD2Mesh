function [vertn, facen] = withTS(vertn, facen)

% translate the model to origin
center = (max(vertn,[],1)+min(vertn,[],1))/2;
vertn = vertn-repmat(center,[size(vertn,1),1]);

% scale the model such that the longest axis is of unit length
bbox = max(vertn,[],1)-min(vertn,[],1);
scale = 1/max(bbox);
vertn = vertn * scale;
end