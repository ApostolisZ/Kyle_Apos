species = {'A', 'B', 'C', 'D'};

Si = [];
Sj = [];
Ss = [];


for irxn = 1:length(1)
rstring = '2 A + B -> C + D';
sep_ind = strfind(rstring,'->');
lhs = strtrim(rstring(1:sep_ind-1));
reactants = strtrim(strsplit(lhs,'+'));
n = regexp(reactants, '(?<stoich>\(?[0-9.]+\)?\s+|)(?<species>\S+)', 'names');
for ireactant = 1:length(n)
    if ~isempty(n{ireactant}.stoich)
        s = regexp(n{ireactant}.stoich, '\(?([0-9.]+)\)?', 'tokens');
        stoich = str2double(s{1});
    else
        stoich = 1;
    end
    ireactant_ind = strfind(species,n{ireactant}.species);
    ireactant_ind = find(~cellfun(@isempty,ireactant_ind));
    Si = [Si;ireactant_ind];
    Sj = [Sj;irxn];
    Ss = [Ss; -stoich];

end


rhs = strtrim(rstring(sep_ind+2:end));
products = strtrim(strsplit(rhs,'+'));
n = regexp(products, '(?<stoich>\(?[0-9.]+\)?\s+|)(?<species>\S+)', 'names');
for iproduct = 1:length(n)
    if ~isempty(n{iproduct}.stoich)
        s = regexp(n{iproduct}.stoich, '\(?([0-9.]+)\)?', 'tokens');
        stoich = str2double(s{1});
    else
        stoich = 1;
    end
    iproduct_ind = strfind(species,n{iproduct}.species);
    iproduct_ind = find(~cellfun(@isempty,iproduct_ind));
    Si = [Si;iproduct_ind];
    Sj = [Sj;irxn];
    Ss = [Ss; stoich];

end




end
S = sparse(Si,Sj,Ss);

foo = 1;

