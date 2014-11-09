function [L_range, y_prime0_n, y_prime0_n1]=SolveBVP
% Here is an Example of solving the Boundary Value Problem by the
% MATLAB build-in function bvp4c or bvp5c    

    %Step1: define the Boundary Value Problem, the ODE and the
    %Boundary condition is defined in the end of this file.
    global A Pr Nb Nt Le M S L    
    S = 1;
    A = 0.5;
    M = 1; % stagnation
    Pr = 0.7;
    Nb = 0.3;
    Nt = 0.7;
    Le = 1.0;
    
    L_init = -1.5; % Take L_init where we have dual-solution.
    L = L_init; 
    guess = make_guess %Step 2: we need to make two initial guess
                       %for L_max
    step_length = 0.1; %you set this to a smaller number.
    guess1 = guess{1}; guess2 = guess{2};
    L_range1 = L_init:-step_length:-2.1; Length1 = length(L_range1);
    L_range2 = L_init:step_length:1;
    L_range = [L_range1, L_range2];
    for j = 1:length(L_range)
        L = L_range( j);
        if j == (Length1 + 1)
            guess1 = guess{1}; guess2 = guess{2};            
        end
        fprintf('%3.6f\n', L);
        %Here we take the solution for the previous L as the
        %initial guess for the next L
        [y1, guess1] = solve_bvp(guess1); 
        [y2, guess2] = solve_bvp(guess2);
        y_prime0_n(j) = y1; 
        y_prime0_n1(j) = y2; 
    end
    [L_range, I] = sort(L_range);
    y_prime0_n = y_prime0_n( I);
    y_prime0_n1 = y_prime0_n1( I);
    plot( L_range, y_prime0_n, 'b-d'); hold on;
    plot( L_range, y_prime0_n1, 'r-*'); hold off;
    
    function guess = make_guess
    % For the given parameter, this function will try to find two
    % different solutions.
        guess1 = ones(7, 1);
        guess2 = ones(7, 1);
        [y1, guess1] = solve_bvp(guess1, 1);
        k = 1;
        guess{1} = guess1; y(k) = y1; k = k + 1;
        for g = 1:2^7
            for i = 1:7
                guess2(i) = 1 - 2*bitand(g, 2^i)/2^i;
            end
            try
                [y2, solu] = solve_bvp(guess2, 1); 
            catch exception
                continue
            end
            count = 0;
            for i = 1:k-1
                if abs( y(i) - y2 ) > 1e-5
                    count = count + 1;
                end
            end
            if count == k - 1
                guess{k} = solu;
                y(k) = y2; k = k + 1
                if k > 2
                    break
                end
            end
        end
    end
    
    function [y, solu ]= solve_bvp(guess, opt)
        infinity = 5; 
        if nargin > 1
            initial_guess = bvpinit( linspace(1e-7, infinity, 100), ...
                                 guess);
        else
            initial_guess = guess;
        end
        options = bvpset('RelTol', 1e-9, 'AbsTol', 1e-9);
        solu = bvp4c(@bvpex, @bcex, initial_guess, options);
        y = solu.y(3, 1);
    end
    
    function ysolu = bvpex(x, y)
    %Here is the ODE.
        ysolu = [y( 2 );
                 y( 3 );
                 - y( 1) * y( 3) + y( 2)*y( 2)  + M * y( 2 )- M * A - A * A;
                 y(5);
                 - Pr * (y(1)*y(5)  + Nb * y(5) * y(7)+ Nt * y(5)*y(5));
                 y(7);
                 - Le * y(1) * y(7) + Nt* Pr * (y(1) * y(5) + Nb * y(5) * y(7)+ Nt * y(5)*y(5))/Nb];
    end
    
    function res = bcex(y0, yinf)
    %Here is the boundary condition of the ODE.
        res = [y0(1) - S;
               y0(2) - L;
               y0(4) - 1;
               y0(6) - 1;
               yinf(2) - A;
               yinf(4);
               yinf(6)];
    end

end
    
