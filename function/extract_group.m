function model = extract_group(GTmodel,group_index)
model = GTmodel;
model.dim(3) = length(group_index);
model.ind_VAR = GTmodel.ind_VAR(group_index);
model.ind = GTmodel.ind(group_index);
model.ind_common = GTmodel.ind_common;
model.ind_differential = GTmodel.ind_differential(group_index);
model.VAR_spectrum = GTmodel.VAR_spectrum(group_index);
model.GC = GTmodel.GC(group_index);
model.A = GTmodel.A(:,:,:,group_index);
model.seed = GTmodel.seed;
end