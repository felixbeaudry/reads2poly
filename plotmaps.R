library(data.table)
library(sp) 
library(maptools) 
library(maps) 
library(mapdata)
library(sfsmisc) 
library(mapproj) 
library(raster)
library(rgeos) 
library(rgdal) 
library(scales)
library(ggmap)
library(mapplots)
library(RgoogleMaps)

gpclibPermit() 

setwd('~/Google Drive/Research/Data/')


pop <- fread("pickupPopsTrim.txt")

USA <- getData("GADM",country="USA",level=1)
states<-c("Texas", "Oklahoma","Louisiana","Alabama","Florida","Georgia","South Carolina","North Carolina","Mississippi","Arkansas")

USA.states <- USA[USA$NAME_1 %in% states,]

state<-geocode(states)
state.x<-state$lon
state.y<-state$lat

USA.bbox <- bbox(USA.states)
xlim <- c(min(USA.bbox[1,1]), max(USA.bbox[1,2]))
ylim <- c(min(USA.bbox[2,1]),max(USA.bbox[2,2]))
plot(USA.states, xlim=xlim, ylim=ylim)


#collection points

add.pie(z=c(1), x=-83.14,y=31.28,  col=c(alpha("black")), labels="", radius=(1/2))
add.pie(z=c(1), x=-81.17,y=31.50,  col=c(alpha("black")), labels="", radius=(1/2))
add.pie(z=c(1), x=-81.5,y=32.27, col=c(alpha("black")), labels="", radius=(1/2))
add.pie(z=c(1), x=-83.4,y=30.34,  col=c(alpha("black")), labels="", radius=(1/4))
add.pie(z=c(1), x=-82.3,y=29.4, col=c(alpha("black")), labels="", radius=(1/4))
add.pie(z=c(1), x=-85.50,y=31.4,  col=c(alpha("black")), labels="", radius=(1/4))
add.pie(z=c(1), x=-86.59,y=31.4, col=c(alpha("black")), labels="", radius=(1/4))
add.pie(z=c(1), x=-81.26,y=34.6, col=c(alpha("black")), labels="", radius=(1/2))
add.pie(z=c(1), x=-79.29,y=34.10, col=c(alpha("black")), labels="", radius=(1/4))
add.pie(z=c(1), x=-80.48,y=33.15, col=c(alpha("black")), labels="", radius=(1/2))
add.pie(z=c(1), x=-78.32,y=34.58, col=c(alpha("black")), labels="", radius=(1/4))
add.pie(z=c(1), x=-76.52,y=35.31,  col=c(alpha("black")), labels="", radius=(1/4))
add.pie(z=c(1), x=-77.36,y=35.15,  col=c(alpha("black")), labels="", radius=(1/4))
add.pie(z=c(1), x=-93.41,y=31.51, col=c(alpha("black")), labels="", radius=(1/4))
add.pie(z=c(1), x=-93.18,y=30.53,  col=c(alpha("black")), labels="", radius=(1/4))
add.pie(z=c(1), x=-94.47,y=30.41, col=c(alpha("black")), labels="", radius=(2/4))
add.pie(z=c(1), x=-95.48,y=32.11,  col=c(alpha("black")), labels="", radius=(2/4))
add.pie(z=c(1), x=-94.59,y=33.10,  col=c(alpha("black")), labels="", radius=(2/4))
add.pie(z=c(1), x=-96.50,y=33.53,  col=c(alpha("black")), labels="", radius=(1/4))
add.pie(z=c(1), x=-95.24,y=34.9,  col=c(alpha("black")), labels="", radius=(2/4))
add.pie(z=c(1), x=-95.37,y=34.53,  col=c(alpha("black")), labels="", radius=(1/4))
add.pie(z=c(1), x=-95.54,y=31.33,  col=c(alpha("black")), labels="", radius=(1/4))


#autosomal
add.pie(z=c(0,1), x=-85.50,y=31.43,  col=c(alpha("orange"),alpha("blue")), labels="", radius=(4/7)) #ALBRU
add.pie(z=c(0,1), x=-86.59,y=31.4,  col=c(alpha("orange"),alpha("blue")), labels="", radius=(2/7)) #ALBRE
add.pie(z=c(0,1), x=-83.4,y=30.34,  col=c(alpha("orange"),alpha("blue")), labels="", radius=(4/7)) #FLJAS
add.pie(z=c(0,1), x=-85.11,y=30.48,  col=c(alpha("orange"),alpha("blue")), labels="", radius=(4/7)) #FLMAR   4
#add.pie(z=c(0,1), x=-82.3,y=29.4, col=c(alpha("orange"),alpha("blue")), labels="", radius=(4/7)) ######FLHAM
add.pie(z=c(0,1), x=-83.14,y=31.28,  col=c(alpha("orange"),alpha("blue")), labels="", radius=(6/7)) #GAGLA
add.pie(z=c(0,1), x=-81.17,y=31.50,  col=c(alpha("orange"),alpha("blue")), labels="", radius=(3/7)) #GABEL
add.pie(z=c(0,1), x=-81.5,y=32.27,  col=c(alpha("orange"),alpha("blue")), labels="", radius=(6/7)) #GASTA
add.pie(z=c(1,0), x=-93.41,y=31.51,  col=c(alpha("orange"),alpha("blue")), labels="", radius=(3/7)) #LABEN
#add.pie(z=c(1,0), x=-93.18,y=30.53,  col=c(alpha("orange"),alpha("blue")), labels="", radius=(1/7)) ####LADER
add.pie(z=c(0,1), x=-78.32,y=34.58,  col=c(alpha("orange"),alpha("blue")), labels="", radius=(4/7)) #NCROS
add.pie(z=c(0,1), x=-76.52,y=35.31,  col=c(alpha("orange"),alpha("blue")), labels="", radius=(3/7)) #NCBAT
#add.pie(z=c(0,1), x=-77.36,y=35.15,  col=c(alpha("orange"),alpha("blue")), labels="", radius=(6/7)) #####NCKIN
add.pie(z=c(0,1), x=-78.46,y=34.38,  col=c(alpha("orange"),alpha("blue")), labels="", radius=(6/7)) #NCELI   
add.pie(z=c(1,0), x=-95.24,y=34.9,  col=c(alpha("orange"),alpha("blue")), labels="", radius=(6/7)) #OKRAT
add.pie(z=c(1,0), x=-95.37,y=34.53,  col=c(alpha("orange"),alpha("blue")), labels="", radius=(6/7)) #OKBAC
#add.pie(z=c(1,0), x=-96.50,y=33.53,  col=c(alpha("orange"),alpha("blue")), labels="", radius=(1/7)) ####OKWIL
add.pie(z=c(0,1), x=-81.26,y=34.6,  col=c(alpha("orange"),alpha("blue")), labels="", radius=(7/7)) #SCPRO
add.pie(z=c(0,1), x=-79.29,y=34.10,  col=c(alpha("orange"),alpha("blue")), labels="", radius=(6/7)) #SCMAR
add.pie(z=c(0,1), x=-80.48,y=33.15,  col=c(alpha("orange"),alpha("blue")), labels="", radius=(5/7)) #SCBRA
add.pie(z=c(1,0), x=-94.47,y=30.41,  col=c(alpha("orange"),alpha("blue")), labels="", radius=(4/7)) #TXLIV
add.pie(z=c(1,0), x=-95.48,y=32.11,  col=c(alpha("orange"),alpha("blue")), labels="", radius=(5/7)) #TXATH
add.pie(z=c(1,0), x=-94.59,y=33.10,  col=c(alpha("orange"),alpha("blue")), labels="", radius=(1/7)) #TXMTP
add.pie(z=c(1,0), x=-95.54,y=31.33,  col=c(alpha("orange"),alpha("blue")), labels="", radius=(2/7)) #TXOAK
add.pie(z=c(1,0), x=-96.51,y=31.7,  col=c(alpha("orange"),alpha("blue")), labels="", radius=(5/7)) #TXROS   

#haplotype

add.pie(z=c(1,0,0,0), x=-83.14,y=31.28,  col=c(alpha("red"),alpha("blue"),alpha("blue"),alpha("orange")), labels="", radius=(1/2))
add.pie(z=c(1,0,0,0), x=-81.17,y=31.50,  col=c(alpha("red"),alpha("blue"),alpha("blue"),alpha("orange")), labels="", radius=(1/2))
add.pie(z=c(0.5,0.5,0,0), x=-81.5,y=32.27,  col=c(alpha("red"),alpha("blue"),alpha("blue"),alpha("orange")), labels="", radius=(1/2))
add.pie(z=c(1,0,0,0), x=-83.4,y=30.34,  col=c(alpha("red"),alpha("blue"),alpha("blue"),alpha("orange")), labels="", radius=(1/4))
add.pie(z=c(1,0,0,0), x=-82.3,y=29.4,  col=c(alpha("red"),alpha("blue"),alpha("blue"),alpha("orange")), labels="", radius=(1/4))
add.pie(z=c(0,1,0,0), x=-85.50,y=31.4,  col=c(alpha("red"),alpha("blue"),alpha("blue"),alpha("orange")), labels="", radius=(1/4))
add.pie(z=c(0,1,0,0), x=-86.59,y=31.4,  col=c(alpha("red"),alpha("blue"),alpha("blue"),alpha("orange")), labels="", radius=(1/4))
add.pie(z=c(0.5,0.5,0,0), x=-81.26,y=34.6,  col=c(alpha("red"),alpha("blue"),alpha("blue"),alpha("orange")), labels="", radius=(1/2))
add.pie(z=c(0,1,0,0), x=-79.29,y=34.10,  col=c(alpha("red"),alpha("blue"),alpha("blue"),alpha("orange")), labels="", radius=(1/4))
add.pie(z=c(0,1,0,0), x=-80.48,y=33.15,  col=c(alpha("red"),alpha("blue"),alpha("blue"),alpha("orange")), labels="", radius=(1/2))
add.pie(z=c(0,1,0,0), x=-78.32,y=34.58,  col=c(alpha("red"),alpha("blue"),alpha("blue"),alpha("orange")), labels="", radius=(1/4))
add.pie(z=c(0,1,0,0), x=-76.52,y=35.31,  col=c(alpha("red"),alpha("blue"),alpha("blue"),alpha("orange")), labels="", radius=(1/4))
add.pie(z=c(0,1,0,0), x=-77.36,y=35.15,  col=c(alpha("red"),alpha("blue"),alpha("blue"),alpha("orange")), labels="", radius=(1/4))
add.pie(z=c(0,0,0,1), x=-93.41,y=31.51,  col=c(alpha("red"),alpha("blue"),alpha("blue"),alpha("orange")), labels="", radius=(1/4))
add.pie(z=c(0,0,0,1), x=-93.18,y=30.53,  col=c(alpha("red"),alpha("blue"),alpha("blue"),alpha("orange")), labels="", radius=(1/4))
add.pie(z=c(0,0,0,1), x=-94.47,y=30.41,  col=c(alpha("red"),alpha("blue"),alpha("blue"),alpha("orange")), labels="", radius=(2/4))
add.pie(z=c(0,0,0,1), x=-95.48,y=32.11,  col=c(alpha("red"),alpha("blue"),alpha("blue"),alpha("orange")), labels="", radius=(2/4))
add.pie(z=c(0,0,0,1), x=-94.59,y=33.10,  col=c(alpha("red"),alpha("blue"),alpha("blue"),alpha("orange")), labels="", radius=(2/4))
add.pie(z=c(0,0,0,1), x=-96.50,y=33.53,  col=c(alpha("red"),alpha("blue"),alpha("blue"),alpha("orange")), labels="", radius=(1/4))
add.pie(z=c(0,0,0,1), x=-95.24,y=34.9,  col=c(alpha("red"),alpha("blue"),alpha("blue"),alpha("orange")), labels="", radius=(2/4))
add.pie(z=c(0,0,0,1), x=-95.37,y=34.53,  col=c(alpha("red"),alpha("blue"),alpha("blue"),alpha("orange")), labels="", radius=(1/4))
add.pie(z=c(0,0,0,1), x=-95.54,y=31.33,  col=c(alpha("red"),alpha("blue"),alpha("blue"),alpha("orange")), labels="", radius=(1/4))


# structure pie charts

add.pie(z=c(0.932,0.011,0,0.057), x=-83.14,y=31.28,  col=c(alpha("red"),alpha("yellow"),alpha("blue"),alpha("black")), labels="", radius=(1/2))
add.pie(z=c(0.998,0.002,0,0), x=-81.17,y=31.50,  col=c(alpha("red"),alpha("yellow"),alpha("blue"),alpha("black")), labels="", radius=(1/2))
add.pie(z=c(0.459,0.469,0,0.071), x=-81.5,y=32.27,  col=c(alpha("red"),alpha("yellow"),alpha("blue"),alpha("black")), labels="", radius=(1/2))
add.pie(z=c(0.950,0,0,0.049), x=-83.4,y=30.34,  col=c(alpha("red"),alpha("yellow"),alpha("blue"),alpha("black")), labels="", radius=(1/4))
add.pie(z=c(0.962,0,0,0.038), x=-82.3,y=29.4,  col=c(alpha("red"),alpha("yellow"),alpha("blue"),alpha("black")), labels="", radius=(1/4))
add.pie(z=c(0,0.831,0,0.169), x=-85.50,y=31.4,  col=c(alpha("red"),alpha("yellow"),alpha("blue"),alpha("black")), labels="", radius=(1/4))
add.pie(z=c(0,0.792,0,0.208), x=-86.59,y=31.4,  col=c(alpha("red"),alpha("yellow"),alpha("blue"),alpha("black")), labels="", radius=(1/4))
add.pie(z=c(0.465,0.483,0,0.052), x=-81.26,y=34.6,  col=c(alpha("red"),alpha("yellow"),alpha("blue"),alpha("black")), labels="", radius=(1/2))
add.pie(z=c(0,0.925,0,0.074), x=-79.29,y=34.10,  col=c(alpha("red"),alpha("yellow"),alpha("blue"),alpha("black")), labels="", radius=(1/4))
add.pie(z=c(0.017,0.964,0,0.019), x=-80.48,y=33.15,  col=c(alpha("red"),alpha("yellow"),alpha("blue"),alpha("black")), labels="", radius=(1/2))
add.pie(z=c(0,0.940,0,0.060), x=-78.32,y=34.58,  col=c(alpha("red"),alpha("yellow"),alpha("blue"),alpha("black")), labels="", radius=(1/4))
add.pie(z=c(0.016,0.919,0,0.064), x=-76.52,y=35.31,  col=c(alpha("red"),alpha("yellow"),alpha("blue"),alpha("black")), labels="", radius=(1/4))
add.pie(z=c(0,0.971,0,0.028), x=-77.36,y=35.15,  col=c(alpha("red"),alpha("yellow"),alpha("blue"),alpha("black")), labels="", radius=(1/4))
add.pie(z=c(0,0,0.823,0.177), x=-93.41,y=31.51,  col=c(alpha("red"),alpha("yellow"),alpha("blue"),alpha("black")), labels="", radius=(1/4))
add.pie(z=c(0,0,0.833,0.167), x=-93.18,y=30.53,  col=c(alpha("red"),alpha("yellow"),alpha("blue"),alpha("black")), labels="", radius=(1/4))
add.pie(z=c(0,0,0.796,0.204), x=-94.47,y=30.41,  col=c(alpha("red"),alpha("yellow"),alpha("blue"),alpha("black")), labels="", radius=(2/4))
add.pie(z=c(0,0,0.851,0.149), x=-95.48,y=32.11,  col=c(alpha("red"),alpha("yellow"),alpha("blue"),alpha("black")), labels="", radius=(2/4))
add.pie(z=c(0,0,0.837,0.162), x=-94.59,y=33.10,  col=c(alpha("red"),alpha("yellow"),alpha("blue"),alpha("black")), labels="", radius=(2/4))
add.pie(z=c(0,0,0.796,0.204), x=-96.50,y=33.53,  col=c(alpha("red"),alpha("yellow"),alpha("blue"),alpha("black")), labels="", radius=(1/4))
add.pie(z=c(0,0,0.855,0.145), x=-95.24,y=34.9,  col=c(alpha("red"),alpha("yellow"),alpha("blue"),alpha("black")), labels="", radius=(2/4))
add.pie(z=c(0,0,0.842,0.158), x=-95.37,y=34.53,  col=c(alpha("red"),alpha("yellow"),alpha("blue"),alpha("black")), labels="", radius=(1/4))
add.pie(z=c(0,0,0.851,0.1490), x=-95.54,y=31.33,  col=c(alpha("red"),alpha("yellow"),alpha("blue"),alpha("black")), labels="", radius=(1/4))




#labels
visited<-c("Houston","Dallas","Austin","Orlando")
visit<-geocode(visited)
visit.x <- visit$lon
visit.y <- visit$lat
text(visit.x,visit.y, labels=visited, cex = .7)

points(pop$Longitude, pop$Latitude, pch=16, cex=1.2)