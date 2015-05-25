function matrix2beTranslated = translateBody(matrix2beTranslated,z)  
[nline ncol] = size(matrix2beTranslated);

    for j = 1:nline
        matrix2beTranslated(j,3)= matrix2beTranslated(j,3)+z ;
    end

