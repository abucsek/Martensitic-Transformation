% Created by Ashley Bucsek, Colorado School of Mines, 2014
% Returns the unit vector of vec_in

function vec_out=unit(vec_in)

if norm(vec_in)==0
    vec_out=vec_in;
else
    vec_out=vec_in/norm(vec_in);
end