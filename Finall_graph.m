close all
clc
clear


%data = cell2mat(struct2cell(load('PLV_alpha')));
data = cell2mat(struct2cell(load('PLV_beta')));
%data = cell2mat(struct2cell(load('PLV_gamma')));
plv_matrix = squeeze(mean(data,1));
%%
%clustering
[acc, c] = weighted_avgClusteringCoefficient(plv_matrix);

%
left = [9 16 17 22:24 27:30 32:79 81:87 89:107 109:111 119 189:192];
right = [1:7 10:14 18:20 25 113:118 121:188];
%
c_right = c(right);
c_left = c(left);
asymmetry = (sum(c_right)-sum(c_left))/(sum(c_right)+sum(c_left))
%%
%global effeicency
data = cell2mat(struct2cell(load('PLV_alpha')));
%data = cell2mat(struct2cell(load('PLV_beta')));
%data = cell2mat(struct2cell(load('PLV_gamma')));
plv_matrix = squeeze(mean(data,1));


n = 192;
sp_plv = sparse(plv_matrix);
G = digraph(plv_matrix);
Gw = G.Edges.Weight;
Dtotal = graphallshortestpaths(sp_plv,'directed',false,'Weights',Gw);
Wn = ((sum(1./(Dtotal+eye(n))))./(n*(n-1)))';
WeightedEG = (sum(sum(1./(Dtotal+eye(n)))) - n)/(n*(n-1))

%
left = [9 16 17 22:24 27:30 32:79 81:87 89:107 109:111 119 189:192];
right = [1:7 10:14 18:20 25 113:118 121:188];
%
Wn_right = Wn(right);
Wn_left = Wn(left);
asymmetry = (sum(Wn_right)-sum(Wn_left))/(sum(Wn_right)+sum(Wn_left))
%%
%centrality

binary_matrix = plv_matrix;
binary_matrix(binary_matrix>=mean(binary_matrix,'all')) = 1;
binary_matrix(binary_matrix<mean(binary_matrix,'all')) = 0;

%%

G = graph(binary_matrix);
n = numnodes(G);

Centrality_degree = centrality(G,'degree')/((n-2)*(n-1)/2);
Centrality_betweenness = centrality(G,'betweenness')/((n-2)*(n-1)/2);

left = [9 16 17 22:24 27:30 32:79 81:87 89:107 109:111 119 189:192];
right = [1:7 10:14 18:20 25 113:118 121:188];

Centrality_degree_right = Centrality_degree(right);
Centrality_degree_left = Centrality_degree(left);
degree_asymmetry = (sum(Centrality_degree_right)-sum(Centrality_degree_left))/(sum(Centrality_degree_right)+sum(Centrality_degree_left));


Centrality_betweenness_right = Centrality_betweenness(right);
Centrality_betweenness_left = Centrality_betweenness(left);
betweenness_asymmetry = (sum(Centrality_betweenness_right)-sum(Centrality_betweenness_left))/(sum(Centrality_betweenness_right)+sum(Centrality_betweenness_left));

%%
%close all
clc
clear

% small world
data = cell2mat(struct2cell(load('PLV_alpha')));
%data = cell2mat(struct2cell(load('PLV_beta')));
%data = cell2mat(struct2cell(load('PLV_gamma')));
plv_matrix = squeeze(mean(data,1));

binary_matrix = plv_matrix;
binary_matrix(binary_matrix>=mean(binary_matrix,'all')) = 1;
binary_matrix(binary_matrix<mean(binary_matrix,'all')) = 0;


[acc, c] = avgClusteringCoefficient(binary_matrix);
G = sparse(binary_matrix);
pl = graphallshortestpaths(G,'directed',false);
pl(pl==Inf) = 0;
G_pl = sum(triu(pl),'all')/192;

Niter = 400;

for i = 1:Niter
[CIJ] = makerandCIJ_und(192,7913);
[acc_r(i), c_r] = avgClusteringCoefficient(CIJ);

Gr = sparse(CIJ);
pl_r = graphallshortestpaths(Gr,'directed',false);
pl_r(pl_r==Inf) = 0;
Gr_pl(i) = sum(triu(pl_r),'all')/192;
end


SW = (acc/mean(acc_r,'all'))/(G_pl/mean(Gr_pl))
%%

%box plot

right_hemisphere_sub = [2,3,4,5,7,8,9,10,11,12,14,16,18,22,24,26];
left_hemisphere_sub = [1,6,13,15,17,19,20,21,23,25];

left = [9 16 17 22:24 27:30 32:79 81:87 89:107 109:111 119 189:192];
right = [1:7 10:14 18:20 25 113:118 121:188];

data = table2array(StrokeCohGraph);

data_RightSub_RightHmph = data(right,right_hemisphere_sub);
data_RightSub_LefttHmph = data(left,right_hemisphere_sub);
data_LeftSub_RightHmph = data(right,left_hemisphere_sub);
data_LeftSub_LeftHmph = data(left,left_hemisphere_sub);

data_RightSub_RightHmph_mean = mean(data_RightSub_RightHmph,1);
data_RightSub_LefttHmph_mean = mean(data_RightSub_LefttHmph,1);
data_LeftSub_RightHmph_mean = mean(data_LeftSub_RightHmph,1);
data_LeftSub_LeftHmph_mean = mean(data_LeftSub_LeftHmph,1);

x = [data_RightSub_RightHmph_mean';data_RightSub_LefttHmph_mean'];
g1 = repmat({'right channels'},16,1);
g2 = repmat({'left channels'},16,1);
g = [g1; g2];
figure
boxplot(x,g)
title('Right hemisphere effected// Lempel-Ziv')

x = [data_LeftSub_RightHmph_mean';data_LeftSub_LeftHmph_mean'];
g3 = repmat({'right channels'},10,1);
g4 = repmat({'left channels'},10,1);
g = [g3; g4];

figure
boxplot(x,g)
title('Left hemisphere effected// Lempel-Ziv')
%%
%box plot

right_hemisphere_sub = [2,3,4,5,7,8,9,10,11,12,14,16,18,22,24,26];
left_hemisphere_sub = [1,6,13,15,17,19,20,21,23,25];

smallworld = table2array(smallworldstroke1);

%smallworld_RightSub = smallworld(right_hemisphere_sub,:);
%smallworld_LeftSub = smallworld(left_hemisphere_sub,:);

x = [smallworld(:,1);smallworld(:,2);smallworld(:,3)];
g1 = repmat({'gamma'},26,1);
g2 = repmat({'beta'},26,1);
g3 = repmat({'alpha'},26,1);
g = [g1;g2;g3];
figure
boxplot(x,g)
title('SmallWorld')
%{
x = [smallworld_LeftSub(:,1);smallworld_LeftSub(:,2);smallworld_LeftSub(:,3)];
g1 = repmat({'gamma'},10,1);
g2 = repmat({'beta'},10,1);
g3 = repmat({'alpha'},10,1);
g = [g1;g2;g3];
figure
boxplot(x,g)
title('Leftt hemisphere effected// SmallWorld')
%}

%%
function [acc, c] = weighted_avgClusteringCoefficient(graph)
deg = (size(graph,1)-1)*ones(size(graph,1),1); %Determine node degrees
cn = diag(graph.^(1/3)*triu(graph.^(1/3))*graph.^(1/3)); %Number of triangles for each node

%The local clustering coefficient of each node
c = zeros(size(deg));
c(deg > 1) = 2 * cn(deg > 1) ./ (deg(deg > 1).*(deg(deg > 1) - 1));

%Average clustering coefficient of the graph
acc = mean(c(deg > 1));
end

function [acc, c] = avgClusteringCoefficient(graph)
deg = sum(graph, 2); %Determine node degrees
cn = diag(graph*triu(graph)*graph); %Number of triangles for each node

%The local clustering coefficient of each node
c = zeros(size(deg));
c(deg > 1) = 2 * cn(deg > 1) ./ (deg(deg > 1).*(deg(deg > 1) - 1)); 

%Average clustering coefficient of the graph
acc = mean(c(deg > 1)); 
end

function [CIJ] = makerandCIJ_und(N,K)
ind = triu(~eye(N));
i = find(ind);
rp = randperm(length(i));
irp = i(rp);

CIJ = zeros(N);
CIJ(irp(1:K)) = 1;
CIJ = CIJ+CIJ';         % symmetrize
end