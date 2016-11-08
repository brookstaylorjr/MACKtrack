function data_out = dropImages (data_in, drop_images)

data_out = data_in;
drops = ismember(data_out.CellData(:,1),drop_images);

data_out.CellData(drops,:) = [];

a = fieldnames(data_out.Measurements);

for i = 1:length(a)
   if length(data_out.Measurements.(a{i}))==length(drops)
       data_out.Measurements.(a{i})(drops,:) = [];
   else
       data_out.Measurements.(a{i})(drop_images,:) = [];
   end
end