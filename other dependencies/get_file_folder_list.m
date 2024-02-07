function [itemPath, itemList] = get_file_folder_list(itemFilter, flag)
%% Extract file or folder list based on the itemFilter
% set flag to 1 to get a file list, set flag to 2 to get a folder list
% itemFilter = '*_*\MU*\WS*\*_l.txt' or '*_*\MU*\WS*'
% Example: 
%   itemPath = get_file_folder_list('*_*\MU*\WS*\*_l.h5', 1)
%   itemPath = get_file_folder_list('*_*\MU*\WS*', 2)
if nargin == 1
   flag = 2;
end
% folderFilter = '*_*\MU*\WS*';
itemList = dir(fullfile(itemFilter));
if flag == 1
    itemList = itemList(~[itemList.isdir]);
else
    itemList = itemList([itemList.isdir]);
end

itemCnt = length(itemList);
itemPath = cell(itemCnt, 1);
for k = 1:length(itemList)
    itemPath{k} = [itemList(k).folder, '\', itemList(k).name];
end
%
end