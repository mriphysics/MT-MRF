function size_var = GetSize(var)
    name = getVarName(var);
    size = whos(name);
    size_var = size.bytes;
end

function out = getVarName(~)
    out = inputname(1);
end