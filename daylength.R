#Forsythe, W. C., Rykiel Jr, E. J., Stahl, R. S., Wu, H. I., & Schoolfield, R. M. (1995). A model comparison for daylength as a function of latitude and day of year. Ecological Modelling, 80(1), 87-95.

daylength=function(lat,day_of_year){
P = asin(.39795*cos(.2163108 + 2*atan(.9671396*tan(.00860*(day_of_year-186)))))

D=24-(24/pi)*acos((sin(0.8333*pi/180)+sin(lat*pi/180)*sin(P))/(cos(lat*pi/180)*cos(P)))
return(D)
}