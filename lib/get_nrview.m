function [azel]=get_nrview(view_choice)
% function [azel,viewname]=get_nrview(view_choice) returns standard azimuth and elevation 
% parameters for 3D views of parameter plots as functions of N
% (downsampling) and R (region size, in downsampled checks)
%
%   input: view_choice
%       if not supplied, or 'show', then azel is structure of available views
%       'ask': asks at keyboard
%       'show': shows values but does not ask
%       {'nrplot','nrline','topdown','matlab'}: 
%             nrplot: default view in ffdm_nrplot
%             nrline: view down a line of constant N*R
%             topdown: top-down view
%             matlab: matlab standard 3D projection
%       if none of the above, then matlab standard is returned
%
%   output: azel, azel(1)=azimuth, azel(2)=elevation, per Matlab's view property
%
%    See also:  FFDM_BTC_SPEC_GEN, FFDM_BTC_PLOT_GEN, GET_NRLIST_AVAIL, FFDM_NRPLOT
%
views=struct();
views.nrplot.azel =[-30.0 15.0];
views.nrplot.desc ='default view in ffdm_nrplot';
views.nrline.azel =[ 40.0 10.0];
views.nrline.desc ='view down a line of constant N*R';
views.topdown.azel=[  0.0 90.0];
views.topdown.desc='top-down view onto (N,R) plane';
views.matlab.azel =[-37.5 30.0];
views.matlab.desc ='matlab standard 3D view';
if nargin==0
    azel=views;
    return
end
view_names=fieldnames(views);
if strcmp(view_choice,'ask') | strcmp(view_choice,'show')
    for iview=1:length(view_names)
        fn=view_names{iview};
        disp(sprintf('%1.0f-> %10s (%40s [az el]=[%6.2f %6.2f])',iview,fn,views.(fn).desc,views.(fn).azel));
    end
    if strcmp(view_choice,'show')
        azel=views;
        return;
    end
    iview=getinp('choice','d',[1 length(view_names)]);
    view_use=view_names{iview};
else
    view_use=view_choice;
end
if strmatch(view_use,view_names,'exact')
    azel=views.(view_use).azel;
else
    azel=views.matlab.azel;
    disp(sprintf('Warning: requested view (%s) not recognized',view_use));
end
return
