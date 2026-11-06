function rc(s,  i,r,c){ r=""; for(i=length(s);i>0;i--){ c=substr(s,i,1); r=r (c=="A"?"T":c=="C"?"G":c=="G"?"C":c=="T"?"A":c) } return r }
