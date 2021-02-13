function B = inf2nan(A)
B = A;
B(abs(A)==inf) = nan;
end