## BIOFILM MODEL COUPLED TO DIFFUSION OF SOLUTES INTO SYSTEM

####### LIBRARIES


library(ReacTran)
library(deSolve)
library(FME)
library(lattice)

communitychars<-matrix(0,nrow=species,ncol=vars)
temp.com<-read.table("communitydata.csv", header=TRUE, sep=",")
communitychars[,1]=temp.com$species
communitychars[,2]=temp.com$abundance
communitychars[,3]=temp.com$LHS
communitychars[,4]=temp.com$mumax
communitychars[,5]=temp.com$ks
communitychars[,6]=temp.com$maint
#####################################################################################################
#                         MODEL STRUCTURE
#
#  IB with each microbe aspecies and a coordinate system spec, east,north
#
#
# grid is the size of the diffusion output (15*15)
# subgrid is wHat wexpect to be able to grow NESTED in the grid (100*100)
#
################################################## DIFFUSION #############################################

########################################################################################################
 #                                       DIFFUSION
##########################################################################################################
# this is based on a Soetart function from ReacTran
# there are many more components that could be defined
# key here are the diffusions rates in x,y and z planes
# boundary conditions in 6 faces of box

# shouLd extend this for other nutrients and modelionic status
#
##
#

diffusion3D <- function(t, Y, par) {

  yy    <- array(dim = c(grid, grid, height), data = Y)  # vector to 3-D array
  dY    <- yy                            

  dY <- dY + tran.3D(C = yy,v.x=3,v.y=3,
      C.x.up = 0, C.y.up = BND, C.z.up = 0,
      C.x.down = 0, C.y.down = 0, C.z.down = 0,
      D.x = Dx, D.y = Dy, D.z = Dz,
      dx = dx, dy = dy, dz = dz, full.check = TRUE)$dC
  return(list(dY))
}
# note to set boundary conditons with flow down x is zero up x is input at 10

## parameters  for diffusion model############################
grid=15
subgrid=100
height=5  
#n=15 # gtid dimensions
dy   <- dx <- dz <- 1   # grid division
Dy   <- Dx <- Dz <- 1.1   # diffusion coeff, X- and Y-direction

y  <- array(dim = c(grid, grid, height), data = 0)     # starting concentrations
#
BND   <- matrix(nrow = grid, ncol = height, 5)      # boundary concentration at upstream point

#


################################ testing model DIFFUSION  ######################
 # steady state this is the way Dana did it
  #   diff2.out<-steady.3D(y, func = diffusion3D, parms = NULL, dimens = c(grid, grid, height),time=0, lrw = 30000000, atol = 1e-10,rtol = 1e-10)
  #   y <- array(dim=c(grid, grid, height), data = diff2.out$y)


#diff2.out<-steady.3D(y, func = diffusion3D, parms = NULL, dimens = c(grid, grid, 5),time=0, lrw = 30000000, atol = 1e-10,rtol = 1e-10)
#     y <- array(dim=c(grid, grid, grid), data = diff2.out$y)


#i=1
#filled.contour(y[, , 1], main=i,color.palette = terrain.colors)

#ybot<-y[,,1]
#####################################################################################








######################################################################################
#                               SOURCE POOL
#
######################################################################################
######################################################################################
#
# draw a randomsample from the abundance species curve of the source to colonise 
#
#
sourcepool<-function(communitychars){
species=nrow(communitychars)
indtot<-sum(communitychars[,2])
spec<-sample(indtot,1)
specnum=0
numspec=communitychars[1,2]
if(spec<=numspec)
   {
   specnum=1
    }
 else
   {
     i=2
     while(i <= species)
        {
        numspec=numspec+communitychars[i,2]
        if(spec<=numspec)
           {
           specnum=i
           i=species
           }
        i=i+1
        }  
   }
#cat("spec sample =",spec," species= ",specnum,"\n")
return(specnum)
}

# FUNCTION TEST
# for(i in 1:100) {sourcepool(communitychars)}

#


#######################################################################################
#######################################################################################
#                              COLONISATION
#
#
#######################################################################################
colonise<-function(allmics, grid, subgrid, invrate){
# where are we full too in allmics
# need to count to end of real entries and then start a new microbe
# or pass value to it from a function that does the counting

tempmic<-matrix(0,nrow=invrate,ncol=8)
# get the colony number for any additions
totmic<-sum(allmics[,1])

# start size for cells
size=3
if(nrow(allmics)==0)
  {
  colony=1
  }
else
  {
  colony=max(allmics[,7])+1
  }
height=1
counter=1
for(i in 1:invrate)
    {
    specnum=sourcepool(communitychars)
    #cat(" spec ",specnum,"\n")
    east<-sample(grid*(subgrid-1),1)
    north<-sample(grid*(subgrid-1),1)
    
   cat(" NEW COLONIST  east ",east," north", north,"\n")
    ## function to check if anything there already
    # need to colonise top of colonies
        
    colon<-0
    
    temp<-subset(allmics,allmics[,2]==east & allmics[,3]==north)   
    if(nrow(temp)==0)
      {
      #
      tempmic[counter,1]=1
      tempmic[counter,2]=east
      tempmic[counter,3]=north
      tempmic[counter,4]=1
      tempmic[counter,5]=specnum                  
      tempmic[counter,6]=size
      tempmic[counter,7]=colony
      colony=colony+1     
     # cat(" is ",counter, "\n")
      totmic=totmic+1
      }
     if(nrow(temp)>0)
      {
      height<-sample(100,1)
      # there is something at base line zero find max for heght and put new cell at height +1
      maxht<-max(temp[,4])
      totmic=totmic+1
      tempmic[counter,1]=1
      tempmic[counter,2]=east
      tempmic[counter,3]=north
      tempmic[counter,4]=maxht+1
      tempmic[counter,5]=specnum                  
      tempmic[counter,6]=size
      tempmic[counter,7]=colony      
      colony=colony+1
      totmic=totmic+1
      }
    counter=counter+1
    }





tempmic<-rbind(as.matrix(allmics),tempmic)
return(tempmic)
}

# TEST FUNCTION this and sourcepool work together
# allmics<-colonise(allmics, grid,subgrid, invrate)
# THIS RETURNS THE COLONISES BIND TO ALLMICS



##########################################################################################
##########################################################################################
##                                         MONOD
#
########################################################################################

monod<-function(thiscellpop, food_avail, communitychars) {


# monod growth component


sizecell<-thiscellpop[6]
species<-thiscellpop[1]
print(species)
mumax<-communitychars[species,4]
halfsatcoeff<-communitychars[species,5]


mu=mumax*(food_avail/(halfsatcoeff+food_avail))

sizecell=sizecell+mu


return(sizecell)
}

###  monod(allmics[1,],4,communitychars)

########################################################################################
########################################################################################
#                     GROWTH AND DIVISION
#
###########################################################################################


# this needs to ascertain total number of cells in part of the grid and then appotion to each 
# the amount that they can get
# suggest that those in layer 1 receive an exponentially smaller component than those on edge
# thi smeans knwoing what layer you are in
# let them 
# way to get number in an area is to sum the subsetted data for a y[,,] grid cell
# WHY NOT SEND OVER A SMALL SECTION AND RUN THE GROWTH ECT IN THAT
# note need to define maxsize soemwhere
# what have we been sent size wise includes all cols sent ?
# but if there is nothing but the one colony then the min and max are one and the same
# so when we have one cell as a colony we should put a buffer of 50 around the whole data
# so mineast is 1 and data start at 50 in east and north

# thiscellpop brings the colony and everything near it over ~ everything is stored elsehwere to get
# a colonised mnmatirx (cant be grown into)
# then subset the colony of interest grow it and send back



growth<-function(thiscellpop,allmicround,ybot, communitychars) {
# note row 1 IS THE COLTOTAKE IN THISCELLPOP
#print("IN GROWTH")
#print(nrow(thiscellpop))
#print(nrow(ybot))

# how much food is there for this colony ?
# note we assume the grid is 100 to the side of the y cell
yrangeE<-floor(thiscellpop[2]/100)+1
yrangeN<-floor(thiscellpop[3]/100)+1


cat("yrangE ",yrangeE," yrangeN ", yrangeN,"\n")
cat("cell is ", thiscellpop, "\n")

food_avail<-ybot[yrangeE,yrangeN]

# now translate into the current_pop grid and store if something present
# the following should be input as variables at the begninning
maxsize=5
mortrate=0.07
efficgrow=0.63
minsize=2

## toggle to stop deads growing or anything else
deadyet=0
## SUBSET THE COLONY TO MODEL FROM THISCELLPOP


## food available per cell in the grid # take off the cells there 


     # now do LHS onthis colony

     
     lhs_order<-sample(4)
     print(lhs_order)
      # growth==1; mortality==2; division==3; 4== maintenace
     for(lhs in 1:4)
        {
        if(lhs_order[lhs]==1 & deadyet==0)
            {
           # print("GROWTH ")
            #   GROWTH CONSUMPTION and METABOLITES
            # now grow using MONOD kinetics 
            
            # for simplicity
            if(food_avail >0)
              {            
              unitgrow<-monod(thiscellpop, food_avail, communitychars)
              thiscellpop[6]<-(unitgrow*efficgrow) *3   ## CHECK
              ybot[yrangeE,yrangeN]=ybot[yrangeE,yrangeN]-unitgrow
              if(ybot[yrangeE,yrangeN]<0)
                 {
                 ybot[yrangeE,yrangeN]=0
                 }

              }
            # otherwise thiscellpop[6] leave as is
            # calculate how much used
            
            
            
            }

         # MORTALITY
       if(lhs_order[lhs]==2 & deadyet==0)
           { 
          # print("MORTALITY ")
           poo<-runif(1)
           if(poo<mortrate)
              {
                       
              east1<-thiscellpop[2]
              north1<-thiscellpop[3]
              height1<-thiscellpop[4]
              cat(" Number ",thiscellpop[7]," deaded at east1 ",east1," north1",north1," height1 ",height1,"\n")
              #current_pop[east1,north1,height1]=0   
               thiscellpop[1:8]=0  
              deadyet=1  # so cannot do any of the other LHS      
              }
            } # mortality end
        # MAINTENANCE
       
       if(lhs_order[lhs]==4 & deadyet==0)
          {
          #print("MAINTENANCE")
          sizenow<-thiscellpop[6]
          specnum<-thiscellpop[1]
          specmain<-communitychars[specnum,6]          
          sizenow<-sizenow-specmain
          if(sizenow<0)
             {
             # kill everything
             print("KILLED in maintenance")
             thiscellpop[1:8]=0
             deadyet=1
             }
           else
             {
             thiscellpop[6]=sizenow
             }
          } # maintenance end


        # DIVISION
       if(lhs_order[lhs]==3 & deadyet==0)
           {
           tempop<-matrix(0,nrow=1,ncol=8)
           #print("DIVISION")
           cat(" size ",thiscellpop[6],"\n")
           print(allmicround)
           if(thiscellpop[6]>maxsize)
               {

                print("IN DIVISION")
                # select from allmicround
                numtochose<-NROW(allmicround)
                posnum<-sample(numtochose,1)
                tempop[1,1]=1
                tempop[1,2]=allmicround[posnum,1]
                tempop[1,3]=allmicround[posnum,2]
                tempop[1,4]=allmicround[posnum,3]
                tempop[1,5]=thiscellpop[5]                    
                tempop[1,6]=minsize
                tempop[1,7]=thiscellpop[7]
                tempop[1,8]=0  

                cat(" selected to move to ",allmicround[posnum,],"\n")
                thiscellpop[6]=thiscellpop[6]-minsize  
                print(" here it is oldpop")            
                print(thiscellpop)
                print("new inidivudal tempop")
                print(tempop) 
                thiscellpop<-rbind(thiscellpop,tempop)
                print("combined")
                print(thiscellpop)
                } 
         

              } # end of division
       } # end of LHS
   
## need to remove post division mortalities if parent cell is killed through failing to maintain or survive mortality

thiscellpop<-subset(thiscellpop,thiscellpop[1]>0)

 
return(list(thiscellpop=thiscellpop,ybot=ybot))

}







### test it

#  growth(thiscellpop,allmicround, ybot,communitychars)




#####################################################################################
######################################################################################
#                               IB MODEL MAIN
#
#
######################################################################################
# allmics has all of the individual bacteria SO need to start a random cell to get it going
#
#
#
#
## the working matrix for each bacteria


IB2<-function(allmics,ybot,communitychars) {


# this will send over one cell at a time and return it so should get return and bind on to what

totalpop<-sum(allmics[,1])
poporder<-sample(totalpop)
print("grow in this order")
print(poporder)
 trip=1
for(count in 1:totalpop)
    {
    print("NEXT CELL")
    coltotake<-poporder[count]
    cat("Cell number ",count, " IT IS COLONY NUMBER ", allmics[coltotake,7],"\n")
    thiscellpop<-allmics[coltotake,]  ## TAKES one cell 
    # ok now search for all within plus one and minus one away and send these over too if EMPTY
    myeast<-allmics[coltotake,2]
    mynorth<-allmics[coltotake,3]
    myht<-allmics[coltotake,4]
    # get those around to which I could divide and send a daughter
    cat("poses",myeast," ",mynorth," ",myht,"\n")
    surround_space <- expand.grid(x=(myeast-1):(myeast+1), y=(mynorth-1):(mynorth+1), z=(myht-1):(myht+1))
    surround_space$address <- paste(surround_space$x, surround_space$y, surround_space$z, sep=".")
    temp<-surround_space[!surround_space$address %in% paste(allmics[,2], allmics[,3], allmics[,4], sep="."),]
    allmicround<-cbind(temp[,1],temp[,2],temp[,3])
    allmicround<-subset(allmicround,allmicround[,1]>0 & allmicround[,2]>0 & allmicround[,3]>0)
    print("CALCULATED ALLMICROUND")
    #note it woudl be more effcient to calculate this only if the cell gets beyond division size from a function inside the growth function
    print(allmicround)
    # this is what is here now send to growth 
   
    if(trip>1)
      {
      print("trip>1")
      temp<-growth(thiscellpop, allmicround,ybot, communitychars)
      newmics<-rbind(newmics,temp$thiscellpop)
      ybot<-temp$ybot
     
      }
    if(trip==1)
      {
       print("trip==1")
       temp<-growth(thiscellpop, allmicround,ybot,comunitychars)
       newmics<-temp$thiscellpop
       ybot<-temp$ybot
       trip=2
      }
    print("************************")
    print("**")
   print(allmics)

    # running colony total for this section (becomes new allmics)
    # now bind to store
    }


# new mics is the latest state of the system after complete run trhough
# at end of run send to and over-write allmics
return(list(newmics=newmics,ybot=ybot))
}






#########################################################################################
####                               FULL COUPLED MODEL
#########################################################################################
##### set up starter



species=4
vars=6
grid=15

## the working matrix for each bacteria
# 1=something, east, north, 

subgrid=100
### what is in allmics ?
### allmics is the individual microbe array.
# column 1 is a column of 1s when an animal is present
# column 2 is east, column 3 is north column 4 is heght  col 5 species number; size is 6,    7 is colony number
# 
# 1=ID; 2=east; 3=north; 4=height; 5=species number; 6= size; 7= colony number

communitychars<-matrix(0,nrow=species,ncol=vars)
temp.com<-read.table("communitydata.csv", header=TRUE, sep=",")
communitychars[,1]=temp.com$species
communitychars[,2]=temp.com$abundance
communitychars[,3]=temp.com$LHS
communitychars[,4]=temp.com$mumax
communitychars[,5]=temp.com$ks
communitychars[,6]=temp.com$maint


##############################STARTER INDIVIDUAL########################################
# THE SPECIES LIFE HISTORY DATA

# note to set boundary conditons with flow down x is zero up x is input at 10

## parameters  for diffusion model############################
grid=15
subgrid=100
height=5  
#n=15 # gtid dimensions
dy   <- dx <- dz <- 1   # grid division
Dy   <- Dx <- Dz <- 1   # diffusion coeff, X- and Y-direction

y  <- array(dim = c(grid, grid, height), data = 0)     # starting concentrations
#
BND   <- matrix(nrow = grid, ncol = height, 2)      # boundary concentration at upstream point

 
# SET UP A RANODM STARTING COORDINATE TO CREAT ALLMICS
### what is in allmics ?
### allmics is the individual microbe array.
# column 1 is a column of 1s when an animal is present
# column 2 is east, column 3 is north column 4 is heght  col 5 species number; size is 6,    7 is colony number




# GET RANODOM POSITIONS AND SPECIES AND SIZE
estart<-sample(grid*(subgrid-1),1)
nstart<-sample(grid*(subgrid-1),1)
#sizestart<-sample(6,1)
sizestart=6
speco<-sample(species,1)
invrate=10
runno=10
allmics<-matrix(data=c( 1, estart, nstart,  1,   speco,   sizestart,   1,   0), nrow=1, ncol=8, byrow=TRUE)

# create storage

results<-matrix(0, nrow=runno,ncol=species+2)
# put food into last and time into first



for(i in 1:runno)
  {
    # DIFFUSSION AND ADVECTION
     cat("model run number  i ****",i,"\n")
     # steady state this is the way Dana did it
     diff2.out<-steady.3D(y, func = diffusion3D, parms = NULL, dimens = c(grid, grid, 5),time=0, lrw = 30000000, atol = 1e-10,rtol = 1e-10)
     y <- array(dim=c(grid, grid, 5), data = diff2.out$y)
     # integration bit more sensitive and can generat goobledegook
     #diff.out <- ode.3D(y, func = diffusion3D, parms = NULL, dimens = c(grid, grid, 5),times = seq(1,2,0.1), lrw = 30000000, atol = 1e-10,rtol = 1e-10)
     #y <- array(dim = c(grid, grid, 5), data = diff.out[nrow(diff.out), -1])
     
     
    # PLOT FOOD AVAILABLE FOLLOWING DIFFUSION
    ybot<-y[,,1]    
    filled.contour(ybot[], main=i,color.palette = terrain.colors)



     # now let the bacteria grow and consume .
    # COLONISE FROM POOL INTO AVAIALBALE SPACE (IN 3D) 
   
     allmics<-colonise(allmics, grid,subgrid, invrate)


    
  
    # ESSENTIAL IB MODEL AND FOOD CONSUMPTION 
    # call growth for individual Life history processes by bacterium
    temp<-IB2(allmics,ybot,communitychars)
    allmics<-temp$newmics
    ybot<-temp$ybot
       cat("\n")
    if(nrow(allmics)>1)
      {
       # COLLATING RESULTS FOR EACH SPECIES AND FOOD
       results[i,1]=i
       results[i,species+2]=sum(ybot)
       results[i,2]=nrow(subset(allmics,allmics[,5]==1))
       results[i,3]=nrow(subset(allmics,allmics[,5]==2))
       results[i,4]=nrow(subset(allmics,allmics[,5]==3))
       results[i,5]=nrow(subset(allmics,allmics[,5]==4))
      }
     else
       {
       print("WENT EXTINCT")
       i=runno
        }

    # PLOT FOOD LEFT  
     cat("Food after growth ")
    filled.contour(log(ybot+1), main=i,color.palette = terrain.colors)
    # put food back
    y[,,1]<-ybot[,]
   

  }



################################ PLOT RESULTS #################################
xyplot(results[,2]+results[,3]+results[,4]+results[,5]~results[,1], ylab="Abundance of 4 species as they colonise", xlab="Time")
### distribution of individuals in colonies

poo<-summary(as.factor(allmics[,7]))

hist(allmics[,7], xlab="Colony number", ylab="Number in colony", breaks=unique(allmics[,7]))

## abundnace of each species of 
truehist(allmics[,5])

####
attach(allmics)
allmics<-data.frame(allmics)
poo<-allmics[order(allmics[,7]),]
detach()







