%must define the name of the struc as fname

fnamestr = getVarName(fname);

for ll=1:length(fname)
eval([fname(ll).Name '=' fnamestr '(ll).Data;']);
end 

function out = getVarName(var)
    out = inputname(1);
end