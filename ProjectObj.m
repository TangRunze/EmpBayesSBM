function val = ProjectObj(x,x_hat)
% Objective functions in projection optimization problem

val = norm(x - x_hat);

end