function [Neq,Req, Seq]=ffdm_btc_s2nr(N,R,S)
% [Neq,Req,Seq]=ffdm_btc_s2nr(N,R,S) converts single values or lists of the parameters 
% N and R to their equivalents that take into account the parameter S
%
% N and R should be one-dimensional and have the same length, S can match
%   or be a scalar
%
% N: downsampling of the original image pixels
% R: number of downsampled blocks in an edge of the region in which btc stats are accumulated
% S: pixel scale at which texture statistics are processed after (N,R) blocking
%
% The equivalences for S=1, referred to the original image, are
% Neq=N*S
% Req=R/S
% Seq=1
%
%   See also:  FFDM_BTC_CALC_GEN, FFDM_BTC_PLOT_GEN, FFDM_BTC_SCAT_GEN.
%
if length(S)==1
    S=repmat(S,length(N),1);
end
Neq=N(:).*S(:);
Req=round(R(:)./S(:));
if any(Req~=R(:)./S(:));
    disp('warning:  some values of S do not divide R');
end
Seq=ones(size(S));
disp(sprintf(' %f2.0 values of N and R converted to equivalents taking S into account'));
disp('     N     R     S   Neq   Req   Seq');
disp([N(:) R(:) S(:) Neq(:) Req(:) Seq(:)])
Neq=reshape(Neq,size(N));
Req=reshape(Req,size(R));
return
