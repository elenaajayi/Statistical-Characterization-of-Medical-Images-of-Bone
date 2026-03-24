function win=mlis_window_setup(map_size,window_type)
% win=mlis_window_setup(map_size,window_type) set up a 2-D window for FFT tapering or other purposes
%
% map_size: size of window
% window_type (compatible with mlis_opts.if_cosbell and mlis_opts.reg_blend_mode)
%   0 or 'hard' or 'cb0': all 1's
%   1 or 'cb1' cosine bell, denominator map_size-1
%   2 or 'cb2' cosine bell, denominator map_size
%        'gau' Gaussian, exp(-(dist)^2/2*sigma^2), where sigma is map_size/4 (so 2 s.d. of the 1-d Gaussian at the edge of the window)
%
% win: window, as array of size map_size x map_size, values in [0 1]
%
%   See also:  MLIS_RUN_SETUP, MLIS_REG_WEIGHTS_SETUP, FFDM_BTC_CALC_GEN, FFDM_BTC_CALC.
%
win=ones(map_size);
switch window_type
    case {0,'cb0','hard'}
    case {1,'cb1'}            
        cosbell=(1-cos(2*pi*[0:map_size-1]/(map_size-1)))/2;
        win=cosbell'*cosbell;
    case {2,'cb2'}
        cosbell=(1-cos(2*pi*[0:map_size-1]/map_size))/2;
        win=cosbell'*cosbell;
    case {'gau'}
        xvals=4*[[0:map_size-1]-(map_size-1)/2]/map_size;
        gau1d=exp(-xvals.^2/2);
        win=gau1d'*gau1d;
    otherwise
    warning('unknown window type');
end
return
