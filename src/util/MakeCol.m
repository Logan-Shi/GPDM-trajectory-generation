%
% y = MakeCol(x)
%
% Turns x into a column vector.
%
% x -- A vector.
% y -- x is x is already a column vector, x' if it is a row vector.
%
function y = MakeCol(x)
  if size(x, 2) > 1
    y = x';
  else
    y = x;
  end;