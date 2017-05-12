function W = generateW(L,per)
% generate the network W
if (nargin < 3)
    name = 'random';
    show_graph = 0;
    quiet = 1;
    arg_opt1 = round(per*(L-1));
end
[~,P] = GGen(L,per,name,show_graph,quiet,arg_opt1);
W = local_degree(P,1);
end