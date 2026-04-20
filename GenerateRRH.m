function RRH = GenerateRRH(N_RRH, Ub, Lb)
    DD_min = 50;
    confict_max = 200; 
    conflict = 0;
    index = 1;
    RRH = ones(N_RRH, 2);
    while index <= N_RRH
        RRH(index, :) = Lb + rand(1,2).*(Ub-Lb);
        flag = 0;
        if any(sqrt(sum((RRH(index, :) - RRH(1 : index - 1, :)) .^ 2, 2)) < DD_min)
            flag = 1;
            conflict = conflict + 1;
            
        end
        
        if flag == 0
            index = index + 1;
            conflict = 0;
        elseif conflict == confict_max
            index = 1;
            conflict = 0;
        end

    end

end