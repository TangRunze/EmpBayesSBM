function val = ProjectObj(x, xHat)
% Objective functions in projection optimization problem

val = norm(x - xHat);

end