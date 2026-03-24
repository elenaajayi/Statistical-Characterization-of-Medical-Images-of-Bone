function [stat_reorg,perm_use]=ffdm_btc_reorg(stat,view_names,stat_name)
% [stat_reorg,perm_use]=ffdm_btc_reorg(stat,view_names,stat_name) is a utility to mirror-flip statistics of btc coordinates of R-facing views 
%
% input:
% stat: statistic values:  dim 1: arbitrary, dim 2: btc coord, dim 3: view number, dim 4: arbitrary
%    order of coordinates on dim3 assumed to conform to order in dict.codel
% view_names: cell array of view names, typically {'LCC','LMLO','RCC','RMLO'}
% stat_name: name of statistic, or empty -- if non-empty, will produce llog message
%
% stat_reorg:  stat, but with dimension 2 reorganized to flip statistics horizontally
% perm_use: permutation used.
%
% stat_reorg(:,:,iview,:)=stat(:,perm_use,iview,:) if view_names{iview} begins with 'R', otherwise unchanged from stat
%
%  See also:  FFDM_BTC_SCAT_GEN, FFDM_BTC_PLOT_GEN, FFDM_BTC_SCAT, FFDM_BTC_CALC_GEN, FFDM_BTC_CALC,
%  BTC_DEFINE, BTC_VFLIP, 
dict=btc_define;
codel=dict.codel;
btc_n=length(codel);
perm_use=[1:btc_n];
for ilet=1:btc_n
    perm_use(ilet)=find(codel==btc_vflip(codel(ilet)));
end
codel_use=codel(perm_use);
%
stat_reorg=stat;
%
for iview=1:length(view_names)
    if (view_names{iview}(1)=='R') %flip only the R-facing views
        stat_reorg(:,:,iview,:)=stat(:,perm_use,iview,:);
        if (~isempty(stat_name))
            disp(sprintf(' view %5s: reorganized coords %s taken from original coords %s in %s',view_names{iview},codel,codel(perm_use),stat_name));
        end
    end
end %iview
return
