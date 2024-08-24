function tf = constraint(x1,x2,x3,x4,z1,z2,z3,z4,dm)

tf1=  sqrt((x1-x2).^2+(z1-z2).^2)>=dm;
tf2=  sqrt((x2-x3).^2+(z2-z3).^2)>=dm;
tf3=  sqrt((x3-x4).^2+(z3-z4).^2)>=dm;
tf4=  sqrt((x1-x3).^2+(z1-z3).^2)>=dm;
tf5=  sqrt((x1-x4).^2+(z1-z4).^2)>=dm;
tf6=  sqrt((x2-x4).^2+(z2-z4).^2)>=dm;

tf = tf1 & tf2 &  tf3 & tf4 &  tf5 & tf6;
end

