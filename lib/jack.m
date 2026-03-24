function [jbias,jdebiased,jvar,jsem]=jack(jnaive,jdrop)
% [jbias,jdebiased,jvar,jsem]=jack(jnaive,jdrop) calculates jackknife statistics from an array
%
%   jnaive: a scalar or array
%   jdrop: a structure array of the parameters in jnaive, each calculated with "drop one"
%    if jnaive is a scalar, then size(jdrop)=[ndrop 1];
%    if jnaive is a column, then jdrop is a 2-d array with size(jdrop)=[length(jnaive),ndrop]
%    if jnaive has size(jnaive,2)>1 or ndims(jnaive)>2, then size(jdrop)=[size(jnaive), ndrop]
%
%   jbias: the values that must be subtracted from jnaive
%   jdebiased: the debiased values of jnaive (jnaive-jbias)
%   jvar: the jackknife variance estimates
%   jsem:  the jackknife standard errors (sqrt(jvar))
%
% See jackbiasvar.doc for documentation, and c30anal_compsrfm for an example of use
%
%   See also:  JACKSA, JACK_DEMO.
%
jbias=[];jdebiased=[];jvar=[];jsem=[];

if prod(size(jnaive))==1
    if ~(ndims(jdrop)==2 & size(jdrop,2)==1)
        error('incompatible dimensions, first input is scalar');
    end
    djack=1;
    repmat_arg=[size(jdrop,1) 1];
elseif ndims(jnaive)==2 & size(jnaive,2)==1
    if ~(size(jdrop,1)==size(jnaive,1))
        error('incompatible dimensions, first input is column vector');
    end
    djack=2;
    repmat_arg=[1 size(jdrop,2)];
else
    djack=ndims(jnaive)+1;
    if ~all(size(jdrop)==[size(jnaive) size(jdrop,djack)])
        error(sprintf('incompatible dimensions, first input is array of dimension %s',ndims(jnaive)));
    end
    repmat_arg=[ones(1,djack-1) size(jdrop,djack)];
end
ndrop=size(jdrop,djack);
rdmean=mean(jdrop,djack);
%do jackknife calcs
jbias=(ndrop-1)*(rdmean-jnaive);
jdebiased=jnaive-jbias;
jvar=(ndrop-1)*mean((jdrop-repmat(rdmean,repmat_arg)).^2,djack);
jsem=jvar.^(0.5);
return
