function X=mat2cell(x)

    X={};
    L=length(x);
    for i=1:L
        X(i)={x(i)} ;   
    end

end