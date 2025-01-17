function errmsg=checkfields(x,allfields)

%CHECKFIELDS check for valid input fields
%   ERRMSG=CHECKFIELDS(X,ALLFIELDS) checks if the struct variable X has 
%   valid fields. ALLFIELDS is a cell variable containing all possible valid
%   fields of X. Each element of the cell is a string. The check is
%   case insensitive.
%
%   ERRMSG=[] if the fields of X are valid, otherwise ERRMSG is a string
%   containing the first wrong field found by CHECKFIEL
%
%   Author: Paolo Serena, 2021
%   University of Parma, Italy

%    This file is part of Optilux, the optical simulator toolbox.
%    Copyright (C) 2021  Paolo Serena, <serena@tlc.unipr.it>
%
%    Optilux is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 3 of the License, or
%    (at your option) any later version.
%
%    Optilux is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

xf = fieldnames(x);
Lx = length(xf);

errmsg=[];
k = 1;
while k <= Lx
    chk = strcmp(xf{k},allfields);
    if ~any(chk)
        errmsg = sprintf('Non existing field %s',xf{k});
        break
    end
    k = k+1;
end
