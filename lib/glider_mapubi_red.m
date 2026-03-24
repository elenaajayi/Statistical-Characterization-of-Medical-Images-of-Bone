function bc_red=glider_mapubired(bc,vecs,vecs_red,g)
% bc_red=glider_mapubi_red(bc,vecs,vecs_red,g) takes the counts of unique block
% configurations on a large glider, and reduces it to counts of unique blocks
% on a subset glider
%
%  bc: block counts on a larger glider, as a row
%  vecs: an array (nvecs x 2) of displacements that define the larger glider
%  vecs_red: an array of displacements that define the smaller glider
%    each row of vecs_red must be a row of vecs
%  g: number of gray levels, defaults to 2
%  
%  bc_red: block counts on the smaller glider
%
%   Conventions for how blocks are packed into bc and bc_red are as in GLIDER_MAPUBI.
%
%  See also:  GLIDER_MAPUBI.
%
if (nargin<=3) g=2; end;
%
r=size(vecs,1);
sumcoords=[1:r]; %which coords to sum over
for irow=1:size(vecs_red,1)
    cmatch=find(min(repmat(vecs_red(irow,:),r,1)==vecs,[],2)==1);
    if length(cmatch)==1
        sumcoords=setdiff(sumcoords,cmatch);
    else
        error('Input arguments vecs and vecs_red are incompatible.');
    end
end
bcp=reshape(bc,repmat(g,1,r));
for coord=sumcoords
    bcp=sum(bcp,coord);
end
bc_red=bcp(:)';
return
