function [mink_wts,mink_off]=get_mink(mink_nns)
% [mink_wts,mink_off]=get_mink(mink_nns) returns the weights and offsets for the Minkowski directions,
% under the convention that 1=black=material, 0=white=empty space, as in Michielsen, K., & De Raedt, H. (2001
% Physics Reports-Review Section of Physics Letters, 347(6), 462–538.
%
% mink_nns: number of nearest neighbors, a vector containing 4, 8 6, -6
%  4: 4-coordinated
%  8: 8-coordinated
%  6: 6-coordinated "L" (connected upper left to lower right)
% -6: 6-coordinated "R" (connected upper right to lower left)
% mink_wts: a matrix of dimension length(nns) x # of btc coords, indicating
%  weighting to make each Minkowski direction (NaN if mink_nns is not in [4 8]
% mink _offs: column of length (nns), offsets
%
% minkowski parameter is mink_off+mink_wts*(btc_coords)'
%
%   See also:  GLM_BTC_MINKDIR, BTC_DEFINE
dict=btc_define;
codel=dict.codel;
btc_n=length(codel);
mink_nn_ref=[4 8 6 -6];
%
%from eqs 1 and 2 of Victor Thengone Rizvi Conte, VR 2015
mink_off_ref=[1 -1 0 0]/16; %amounts to add to Minkowski weights*btc coords to get Minkowski params
mink_wts_ref(1,:)=[-4 -2 -2 +1 +1 +1 +1 +1 +1 +1]/16; %coefs for g,b,c,d,e,t,u,v,w,a
mink_wts_ref(2,:)=[-4 +2 +2 -1 -1 +1 +1 +1 +1 -1]/16; %coefs for g,b,c,d,e,t,u,v,w,a
mink_wts_ref(3,:)=[-2  0  0  0  0  0 +1  0 +1  0]/8; %nonzero coefs for g, u, w
mink_wts_ref(4,:)=[-2  0  0  0  0  +1 0 +1  0  0]/8; %nonzero coefs for g, t, v
%
nmink=length(mink_nns);
mink_wts=nan(nmink,btc_n);
mink_off=nan(nmink,1);
for imink=1:nmink
    p=find(mink_nns(imink)==mink_nn_ref);
    if ~isempty(p)
        mink_wts(imink,:)=mink_wts_ref(p,:);
        mink_off(imink,:)=mink_off_ref(p);
    end
end
return
