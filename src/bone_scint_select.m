function [xlsfile,desc,opts]=bone_scint_select(setups)
%[xlsfile,desc,opts]=bone_scint_select(setups) chooses a database of images and choices for how to read it
%
% setups: optional input with nonstandard setups
%
% xlsfile: default name of xls file with images
% desc: text descriptor of database and options
% opts: optoins suitable for bone_xls_read
%
%   See also: BONE_XLS_READ, BONE_BTC_DEMO, BONE_DBASE_DEMO, BONE_PSPEC_DEMO.
%
if (nargin<1)
    setups=cell(0);
    setups{1}.desc='bone radiographs from MedPix database';
    setups{1}.xlsfile='bonedatabase.xlsx';
    setups{1}.nroi_col_tag='how many ROI?';
    setups{1}.roi_infix='';
    %
    setups{2}.desc='bone radiographs from MedPix database, interiors only';
    setups{2}.xlsfile='bonedatabasesubroi.xlsx';
    setups{2}.nroi_col_tag='how many iroi';
    setups{2}.roi_infix='i';
    %
    setups{3}.desc='scintigraphy from MedPix database';
    setups{3}.xlsfile='scindatabase.xlsx';
    setups{3}.nroi_col_tag='how many ROI?';
    setups{3}.roi_infix='';
end
for isetup=1:length(setups)
    disp(sprintf('setup %1.0f: %s, file name %s',isetup,setups{isetup}.desc,setups{isetup}.xlsfile));
end
isetup=getinp('choice','d',[1 length(setups)]);
xlsfile=setups{isetup}.xlsfile;
desc=setups{isetup}.desc;
opts=rmfield(setups{isetup},{'xlsfile','desc'});
return

    