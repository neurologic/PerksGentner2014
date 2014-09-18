function norm_data = normdata(data)

data = data - min(data);
norm_data = data ./ max(data);