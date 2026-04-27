function [si_X,true_labs] = info_data(num)

if num == 1

    load("Data_Buettner.mat");

elseif num == 2

    load("Data_Deng.mat");
    
elseif num == 3

    load("Data_Schlitzer.mat");

elseif num == 4

    load("Data_Macosko.mat");
    
    J = find(true_labs<=4);
    
    true_labs = true_labs(J);
    
    in_X = in_X(J,:);
        
elseif num == 5

    load("Data_Pollen.mat");

elseif num == 6

    load("Data_Tasic.mat");
    
    in_X = double(in_X);
    
    J = find(true_labs<=9);
    
    true_labs = true_labs(J);
    
    in_X = in_X(J,:);
    
elseif num == 7

    load("Data_Ting.mat");

elseif num == 8

    load("Data_Treutlin.mat");

elseif num == 9

    load("Data_Zeisel.mat");

end

si_X = func_simlarity(in_X);

end