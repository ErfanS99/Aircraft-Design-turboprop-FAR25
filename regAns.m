function W_E1 = regAns(W_TO_guess)
    A = 0.108637;
    B = 1.050745;
    W_E1 = 10 .^ ((log10(W_TO_guess) - A) / B);
end
