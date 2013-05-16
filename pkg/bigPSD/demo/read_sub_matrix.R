n=9229
p=500568
from_column=1
to_column=10000
centers=rep(0,p)
weights=rep(1,p)
file="/Users/paulino/Desktop/genosInt.bin"

a=read_sub_matrix(n,p,from_column,to_column,file,centers,weights)
