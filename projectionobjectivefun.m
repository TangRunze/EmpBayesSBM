function val = projectionobjectivefun(x, xHat)
% Objective functions in projection optimization problem

val = norm(x - xHat);

end