function [jbias,jdebiased,jvar,jsem]=jacksa(jnaive,jdrop)
% [jbias,jdebiased,jvar,jsem]=jacksa(jnaive,jdrop) calculates jackknife statistics from the fields 
% in a structure
%
%   jnaive: a structure with one or more fields; each field can be an array
%   jdrop: a structure array of the parameters in jnaive, each calculated with "drop one"
%     All fields in jnaive must be present in jdrop
%
%   jbias: the values that must be subtracted from jnaive
%   jdebiased: the debiased values of jnaive (jnaive-jbias)
%   jvar: the jackknife variance estimates
%   jsem:  the jackknife standard errors (sqrt(jvar))
%
% See jackbiasvar.doc for documentation, and c30anal_compsrfm for an example of use
%
%   See also:  BOOTSA, JACK.
%
jbias=[];jdebiased=[];jvar=[];jsem=[];
ndrop=length(jdrop);
names=fieldnames(jnaive);
for iname=1:size(names,1);
   u=deblank(names(iname,:));
   fname=u{1}; %convert from cell array to string
   val=getfield(jnaive,fname);
   shape=size(val);
   r=reshape(val,1,prod(shape));
   clear rd
   for k=1:ndrop
      val=getfield(jdrop(k),fname);
      rd(k,:)=reshape(val,1,prod(shape));
   end
   rdmean=mean(rd,1);
   %do jackknife calcs
   bias=(ndrop-1)*(rdmean-r);
   debiased=r-bias;
   var=(ndrop-1)*mean((rd-repmat(rdmean,ndrop,1)).^2,1);
   sem=var.^(0.5);
   jbias=setfield(jbias,fname,reshape(bias,shape));
   jdebiased=setfield(jdebiased,fname,reshape(debiased,shape));
   jvar=setfield(jvar,fname,reshape(var,shape));
   jsem=setfield(jsem,fname,reshape(sem,shape));
end
return
