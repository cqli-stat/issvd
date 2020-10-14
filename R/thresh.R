### z: argument
### type: thresholding rule type
### delta: thresholding level
### a: default choice for SCAD penalty

thresh <- function(z, type, delta, a=3.7){

  if(type=="soft"){
    return(sign(z)*(abs(z)>=delta)*(abs(z)-delta))
  }
  if(type=="hard"){
    return(z*(abs(z)>delta))
  }
  if(type=="SCAD"){
    return(sign(z)*(abs(z)>=delta)*(abs(z)-delta)*
             (abs(z)<=2*delta)+((a-1)*z-sign(z)*a*delta)/(a-2)*
             (2*delta<abs(z))*(abs(z)<=a*delta)+z*(abs(z)>a*delta))
  }
}
