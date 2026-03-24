function cns=glider_lkpind(inds,gs)
% cns=glider_lkpind(inds,gs) looks up one or more indices, to determine their check
% number in the glider
%
% inds: an array with 2 columns, x-index and y-index (1-based)
% gs:  a glider structure
% gs.matrix: array of 0's and 1's
% gs.size:  size of bounding box (rows and cols)
% gs.rank:  number of occupied cells
% gs.inds:  pointers (1-based) into occupied cells; size(g.inds)=[g.rank 2]
%
% cns: a column of length equal to size(inds,1), indicating which check
% number (1 to gs.rank) corresponds to each check
%
% if the check is not found, cns=0
%
indvec=(inds(:,1)-1)*gs.size(2)+inds(:,2);
gsvec=(gs.inds(:,1)-1)*gs.size(2)+gs.inds(:,2);
for k=1:length(indvec)
    if (ismember(indvec(k),gsvec))
        cns(k,1)=min(find(indvec(k)==gsvec));
    else
        cns(k,1)=NaN;
    end
end
return

