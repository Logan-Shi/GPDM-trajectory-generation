%
% y = MakeRow(x)
%
% Turns x into a row vector.
%
% x -- A vector.
% y -- x is x is already a row vector, x' if it is a column vector.
%
function y = MakeRow(x)
  if size(x, 1) > 1
    y = x';
  else
    y = x;
  end;
