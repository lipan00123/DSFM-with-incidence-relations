clear;
load 020_smallplant_wts
load 020_smallplant_p500
mex RCDM_cversion.cpp;
mex RCDM_cversion_greedypar.cpp;
mex ACDM_cversion.cpp;
mex APcompact_cversion.cpp;
mex APfull_cversion.cpp;
%% preprocessing/preparing data
%unary potential
a = reshape(unary,[427,640]);
%number of columns
col = 640;
%number of rows
row = 427;
%number of superpixels
region = 500;
%number of pixels
N = length(unary);
%number of decomposed parts of the whole submodular function
R = col-1 + row-1 + region;
%incidence_list{i} is the incidence set of the ith submodular function 
incidence_list = cell(R,1);
%parameter_list{i} is the parameter for the ith submodular function 
parameter_list = cell(R,1);
%submodular_type{i} is the type the ith submodular function 
submodular_type = cell(R,1);
%Submodular functions in columns can be decomposed into parallel edges.
%submodular_type is "para_edge"
for i = 1:col-1,
    for j = 1:row,
        incidence_list{i} = [incidence_list{i} POS(j,i,row) POS(j,i+1,row)];
        parameter_list{i} = [parameter_list{i} 0.08 * W2(j,i)];
        submodular_type{i} = 'para_edge';
    end
end
%Submodular functions in rows can be decomposed into parallel edges.
%submodular_type is "para_edge"
for i = 1:row-1,
    for j = 1:col,
        incidence_list{i+col-1} = [incidence_list{i+col-1} POS(i,j,row) POS(i+1,j,row)];
        parameter_list{i+col-1} = [parameter_list{i+col-1} 0.08 * W1(i,j)];
        submodular_type{i+col-1} = 'para_edge';
    end
end
%Submodular functions in superpixels can be decomposed into F(S) = |S||e/S|.
%submodular_type is "concave_card", parameter_list{i}[j] = F([j]) - F([j-1])
for i = 1:region,
    subset = find(sup_img == i);
    card = length(subset);
    subset_row = rem(subset-1, row) + 1;
    subset_col = (subset-1 - rem(subset-1, row))/row + 1;
    incidence_list{i+col-1+row-1} = reshape(POS(subset_row, subset_col, row),1, card);
    parameter_list{i+col-1+row-1} = 0.008*(card + 1 - 2*(1:card));
    submodular_type{i+col-1+row-1} = 'concave_card';
end
bias_vec = reshape(a,1,N);

% K = round(alpha*R) many decomposed parts are simulated to be computed parallelly 
alpha = 0.02;
K = round(alpha*R);
% T is the number iterations of parallel RCD, ACD
T = 600/alpha;
% record_dis controls how many iterations to perform between recording two
% gaps
record_dis = T/200;
% c is an empirical parameter for ACD to control the number of iterations
% for each outloop
c = 10;

%to obtain greedy partition
partitionlabel = [];
subdegree = inf;
for i = 1:10,
    [temppartition, tempsumdegree] = greedy_partition(incidence_list, N, R, K);
    if tempsumdegree < subdegree,
        subdegree = tempsumdegree;
        partitionlabel = temppartition;
    end
end
G_part = max(partitionlabel);
for i = 1:G_part,
    partition{i} = (find(partitionlabel == i))';
end

%% solve DSFM problems 
% record contains the continuous duality gaps; 
% record_discrete contains the discrete duality gaps;
%parallel RCD for DSFM
[y, record{1}, record_discrete{1}] = RCDM_cversion(incidence_list, parameter_list, submodular_type, bias_vec, N, R, K, T, record_dis);
%parallel RCD + greedy-partition for DSFM
[y, record{2}, record_discrete{2}] = RCDM_cversion_greedypar(incidence_list, parameter_list, submodular_type, bias_vec, N, R, K, T, record_dis, G_part, partition);
%parallel ACD for DSFM
[y, record{3}, record_discrete{3}] = ACDM_cversion(incidence_list, parameter_list, submodular_type, bias_vec, N, R, K, T, record_dis, c);

% T is the number iterations of parallel AP
T = 600;
record_dis = T/200;

%IAP for DSFM
[y, record{4}, record_discrete{4}] = APcompact_cversion(incidence_list, parameter_list, submodular_type, bias_vec, N, R, T, record_dis);
%AP for DSFM
[y, record{5}, record_discrete{5}] = APfull_cversion(incidence_list, parameter_list, submodular_type, bias_vec, N, R, T, record_dis);


