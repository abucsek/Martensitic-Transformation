function [vec_out1,vec_out2,vec_out3]=reciprocal(vec_in1,vec_in2,vec_in3)

V=dot(vec_in1,cross(vec_in2,vec_in3));

vec_out1=(cross(vec_in2,vec_in3)/V);  
vec_out2=(cross(vec_in3,vec_in1)/V);
vec_out3=(cross(vec_in1,vec_in2)/V);