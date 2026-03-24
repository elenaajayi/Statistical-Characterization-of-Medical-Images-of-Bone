function [trimrange,msg]=roi_trim(mask,downsample)
% [trimrange,msg]=roi_trim(mask,downsample) trims a mask file and determines that the starting phase is consistent 
% with a downsampling scale
% 
% mask: an array of 0's and 1's, 1's indicate the region of interest
% downsample: an integer, typically power of 2
%
% trimrange=[xlo xho ylo yhi], range within mask that contains the entire roi, 
%   with every transition from 0 to 1 occurring at a pixel that is =1 (mod downsample), and every 
%   transition from 1 to 0 occurring at a pixel that is =0 (mod downsample).
% msg: message if mod condition does not hold
%
% If this is not possible then trimrange returned empty
%
%  See also:  BTCSTATS_ROI_TIFF_DEMO.
%
trimrange=[];
msg=[];
for icoord=1:2
    used=find(any(mask,3-icoord));
    trimrange(2*(icoord-1)+[1 2])=[min(used) max(used)];
end
%find the transitions - first transition must occur at low end of range
mask_trimmed=mask(trimrange(1):trimrange(2),trimrange(3):trimrange(4));
for icoord=1:2
    uvec=([1 2]==icoord);
    mask_shifted=mask_trimmed(1+uvec(1):end,1+uvec(2):end);
    mask_noedge=mask_trimmed(1:end-uvec(1),1:end-uvec(2));
    transition10=find(any(mask_shifted<mask_noedge,3-icoord));
    transition01=find(any(mask_shifted>mask_noedge,3-icoord));
    %
    nfail=sum((mod(transition10,downsample)~=0));
    if nfail>0
        msg=strvcat(msg,sprintf('dim %1.0f: boundary condition fails for 1->0 transitions in %2.0f locations',icoord,nfail));
    end
    nfail=sum((mod(transition01,downsample)~=0));
    if nfail>0
        msg=strvcat(msg,sprintf('dim %1.0f: boundary condition fails for 0->1 transitions in %2.0f locations',icoord,nfail));
    end
end
if ~isempty(msg)
    trimrange=[];
end
return

