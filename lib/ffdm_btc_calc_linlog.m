function [patches_xform,label]=ffdm_btc_calc_linlog(patches,ds)
% [patches_xform,label]=ffdm_btc_calc_linlog(patches,ds) is a utility for
% ffdm_btc_calc_gen that determines whether to transform image patches linearly or
% logarithmically, and then carries this out
%
% input:
% patches: image patches, as a multiimensional array
% ds: data source descriptor, with fields roi_range_init and roi_range, giving low and high range
%    before and after transformation
%
% output:
% patches_xform: transformed patches
% label: 'log or 'linear'
%
%   See also:  FFDM_BTC_CALC_GEN.
%
iflog=getinp('1 to transform logarithmically','d',[0 1]);
lum_max=ds.roi_range_init(2); %maximum lum value
patches_xform=min(max(ds.roi_range(1),patches),lum_max);
if (iflog)
    label='log';
    patches_xform=log(patches_xform+1)/log(lum_max+1);
else
    label='linear';
    patches_xform=patches_xform/lum_max;
end
return
