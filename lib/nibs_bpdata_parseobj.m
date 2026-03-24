function obj_table=nibs_bpdata_parseobj(obj_list)
% obj_table=nibs_bpdata_parseobj(obj_list) parses a cell array of object
% descriptors into a struture
%
% obj_list: a cell array with elements such as {'a67-1','a67-2','m69-1'}
% obj_table: a structure with fields 
%    obj_table.type (typically 'a' or 'm')
%    obj_table.pid_num (PID number, 67 or 69 in the above example)
%    obj_table.obj_num (object number, 1 or 2 in the above example)
nobj=length(obj_list);
obj_table=struct;
obj_table.type=repmat(' ',nobj,1);
obj_table.pid_num=zeros(nobj,1);
obj_table.obj_num=zeros(nobj,1);
for iobj=1:nobj
    desc=obj_list{iobj};
    obj_table.type(iobj)=desc(1);
    dashloc=find(desc=='-');
    if length(dashloc)==1
        obj_table.pid_num(iobj,1)=str2num(desc(2:dashloc-1));
        obj_table.obj_num(iobj,1)=str2num(desc(dashloc+1:end));
    end
end
return

