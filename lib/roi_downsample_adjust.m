function roi_data_adj=roi_downsample_adjust(roi_data,downsample)
% roi_data_adj=roi_downsample_adj(roi_data,downsample) adjusts a region-of-interest mask
% so that its blocks align with squares of size downsample
%
% roi_data: the original roi mask, 1 within the region of interest, 0 outside
% downsample: amount of downsampling, typically power of 2
%
% roi_data_adj: roi mask adjusted so that blocks align with downsamping
%
%   See also:  ROI_TRIM, BTC_ROI_TIFF_DEMO
%
%set a block to 1 only if all pixels are 1
%
roi_data_blocked=reshape(roi_data,[downsample,size(roi_data,1)/downsample,downsample,size(roi_data,2)/downsample]);
roi_data_min=squeeze(min(min(roi_data_blocked,[],3),[],1));
roi_data_adj=reppxl(roi_data_min,downsample);
return
end
