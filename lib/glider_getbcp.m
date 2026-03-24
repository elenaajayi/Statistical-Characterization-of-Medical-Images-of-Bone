function [bcr,bca]=glider_getbcp(gs,ng,opts)
% [bcr,bca]=glider_getbcp(gs,ng,opts) gets independent parameters for block configurations
% for a glider specified by gs, and ng gray levels
%
% gs.matrix: array of 0's and 1's
% gs.size:  size of bounding box (rows and cols)
% gs.rank:  number of occupied cells
% gs.inds:  pointers (1-based) into occupied cells; size(g.inds)=[g.rank 2]
%
% ng: number of gray levels, defaults to 0
% opts: options
% opts.getbcp_fill: value to fill irrelevant entries in bcr, bca with, defaults to NaN
% opts.getbcp_show: 1 to show progress, 0 (default) otherwise
%
% bcr:  array of gs.rank columns, one row for each independent parameter
%       for the independent probabilites on subsets of a glider
%       assuming that subsets that differ by a translation have the same
%       probabilities
%       first column: number of occupied cells
%       second column: which subglider this is
%       columns 3 to gs.rank+2: contents of the cell (1 to ng)
%       Note that probabilities involving a contents of 0 can be inferred
%       from above, by consistency conditions with lower ranks
% bca:  similar to bcr, but subsets that differ by a translation may
%       have different probabilities
%       first column: number of occupied cells
%       second column: which subglider this is
%       columns 3 to gs.rank+2: contents of the cell
%       last column: pointer to corresponding row of bcr
%
%   See also:  INT2NARY, GLIDER_ADDSUBS.
%
if (nargin<=1) ng=2; end
if (nargin<=2) opts=[]; end
opts=filldefault(opts,'getbcp_fill',NaN);
opts=filldefault(opts,'getbcp_show',0);
%
bcr=zeros(0,gs.rank+2);
bca=zeros(0,gs.rank+3);
for r=1:gs.rank
    subs=gs.subs{r};
    nsubs=size(subs.relvecs,3);
    if (opts.getbcp_show)
        disp(sprintf(' rank %2.0f has %5.0f subsets',r,nsubs))
    end
    for s=1:nsubs
        offsets=subs.offsets{s};
        for t=1:size(offsets,1)
            coords=subs.relvecs(:,:,s)+repmat(offsets(t,:),r,1);
            cns=glider_lkpind(coords,gs)';
            row=ismember([1:gs.rank],cns);
            if (ng==2)
                if (t==1)
                    bcrptr=size(bcr,1)+1;
                    bcr=[bcr;[r s row]];
                end
                bca=[bca;[r s row bcrptr]];
            else
                nfills=(ng-1)^r; %number of fills of r blocks with [1:ng-1]
                combs=1+int2nary([0:(nfills-1)]',ng-1,r); %a way to make combinations of [1:ng-1]
                rows=zeros(nfills,gs.rank);
                for br=1:r
                    rows(:,cns(br))=combs(:,br);
                end
                if (t==1)
                    bcrptr=size(bcr,1)+[1:nfills]';
                    bcr=[bcr;[repmat([r s],nfills,1) rows]];
                end
                bca=[bca;[repmat([r s],nfills,1) rows bcrptr]];
            end
        end %t
    end %s
end %rank
bcr(find(bcr==0))=opts.getbcp_fill;
bca(find(bca==0))=opts.getbcp_fill;
return


