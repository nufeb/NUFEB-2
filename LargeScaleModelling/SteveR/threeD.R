# IB and PDE mix 

library(ReacTran)
library(deSolve)
library(FME)
library(lattice)


#  Couple a 3 d diffusion model to an IB film over which the 3d water columnis flowing
#
#
# steady state or numerical integrationthrough time cna be done here
#check which to get results out correctly 
#
#







################################################################################################
                               IB MODEL
################################################################################################
# 
# This is a simple IB model that allows bacteria to grow 
# cells occupy a surface in a 3d box of size grid
# the surface is a lattice of size grid. Each lattice positon can support a population (here 100)
#  cells grow from a minimum viable size up to a maximum size before they divide into two
#  
#
# SO
# if population reaches 100 in a cell then they overflow into adjacent cells (queens move)
# note position 101 is the count of bacteria in that lattice position
# note position 102 is the food brought by diffusion to that lattice position
#  note position 103 is the movers leaving that cell to wander lonely
# mortality is imposed at a fixed rate to all bacterial cells
# if the cells drop below a minimum size (maintenance) they also die
#
# note: could use position 1:100 to create layers and then evaluate growth and uptake in each layer
# note: the current format cannot be used with multiple species without extending sizenum array into 4th dimension
# could have microbes as objects with position, size, age and type as attributes
# there will come a point where number of objects will be less efficient to model than individuals within an areal array
#
#


siz_ind<-function(sizenum,grid,row,col,grate,mortality,maxcell,mainenance){
         # this will grow cells and let them split
         # for now embed the growth and survival parameters in the function
         #
         # note hard coded the sizenum matrix to max population (carrying capacity) of 100
         #
         # IB parameters

babycell=maxcell/2

# grate=4  # weight added from resource per time step in full food conditions
 # a killing ofr the newbies ~
 
# it will add to number until a size is reached then one will split
# need to consider doingthis at random through the 100 cells to be true IB OR randomise in replacing back into array after mortality
# 
totcell=0
movers=0
### work out how many cells are there and how much food has been provided by diffusion
oldcell<-sizenum[row,col,101]
avail<-sizenum[row,col,102]
maxusable=grate*oldcell
amount_each<-avail/oldcell   ###<------ can get divide by zero

tempgrate=0 # realised grate given limited food set to zero to start
finavail=0 # amount unused after division and growth returned to next diffusion process set to zero first

# growth up to division  grate is accumulation
# can grow to maxcell in one time step and divide
# BUT only one divide. If the division occurs sizes set to babycell and balance of adult minus baby 
# need to dcrement availability by the amount needed to get to divide point  
cat("amount each ",amount_each,"  cells to eat it ",oldcell,"\n")
## randomly kill some 
# inefficient looping but need to keep lists of cells with no gaps, do some maintenance, too
# temporary storage matrx to allow for random mortality
temparr<-matrix(0,100,1)
counter=0
# count through all of the cells and draw a random number to make a decision on killing them
for(i in 1:oldcell)
   {
   killer=runif(1)
   #cat("killer ",killer,"\n")
   if(killer<=mortality)
     {
     #cat("I survived \n")
     # it lives for now
     # doesit have enough to maintain itself remove maintenance for tim step and kill if necessary
     sizenum[row,col,i]=sizenum[row,col,i]-maintenance
    # cat("doing some death to number ",i,"\n")
     if(sizenum[row,col,i]>=maintenance)
       {
        counter=counter+1
        temparr[counter,1]=sizenum[row,col,i]   
        # cat("I lived ",i,"\n")
        } 
      else
        {
        # cat(" I died ",i,"\n")
        }

      } 

   }
# this has left holes in the sizenum matrix where cells have died therefore re-sort
cat("after mort total alive ",counter," left from ", oldcell,"\n")
sizenum[row,col,103]=(oldcell-counter)
# clear population array at this grid cell
sizenum[row,col,]=0
# put back cells ie shuffle verything back to start point in matrix to avoid having zeroes in mid column
# should do this at random 
for(i in 1:counter)
  {
  sizenum[row,col,i]=temparr[i,1]
  }
cat("old ",oldcell," new",counter,"\n")
oldcell=counter  # because we killed some of oldcell by the mortality factor 
# only do some growth and movement if there is something to grow after the mortality
if(oldcell>0)
  {
  for(i in 1:oldcell)
    {
    tempgrate=0
    ### look at how much there is for the population
    if(amount_each<=grate)                                      
       {
       tempgrate=amount_each
       finavail=0
       }
    if(amount_each>grate)
       {
       tempgrate=grate
       finavail=avail-(oldcell*grate)
       }
    # save food left if any
    sizenum[row,col,102]=finavail
    # increment cell size by estimated growth rate
    sizenum[row,col,i]=sizenum[row,col,i]+tempgrate
    #cat(" * ", sizenum[row,col,i],"\n")
    # if cell bigger than macell size divide
    if(sizenum[row,col,i]>maxcell)
       {
      sizenum[row,col,i]=sizenum[row,col,i]-babycell # note this leaves parent slightly bigger which is realistic 
       # cat(" * ", sizenum[row,col,i],"\n")
       totcell=totcell+1
      
       }
     }
 


   # now calculate how much space there is at this lattice point
   # increment the population withthe babies and store their size

  cat(" cells reprdoduced ",totcell,"\n")
  ## add new cells up to 100 (carrying capacity)
 if(oldcell+totcell<=100)
    {
    for(i in 1:totcell)
       {
       sizenum[row,col,oldcell+i]=babycell
       }
    sizenum[row,col,101]=oldcell+totcell
    sizenum[row,col,102]=0
    }  
  cat(" col together ",(sizenum[row,col,101]+totcell),"\n")
   supertot<-(sizenum[row,col,101]+totcell)
   # if we have more cells than 100 we must move some 
   if(supertot>100)
     {
     # we have some surplus population as the population can only go to 100
     for(q in sizenum[row,col,101]:100)
        {
        #cat(" q ",q," ol ",oldcell,"\n")
        sizenum[row,col,q]=babycell
      
        }
    movers<-supertot-100
    sizenum[row,col,101]=100
   sizenum[row,col,102]=avail
   # number that have to move
   
   cat(" movers ",movers,"\n")
   }    
   cat(" row ",row," col ",col,"\n")
   ## to move them they can go +1 and -1 in row and col select where moved to at random
   ## note: should turn what follows into a dispersal function
  if(movers>0)
    {
   # must check where anywhere is feasible locally
   #  there is an initial check of whether or not space is available using fullup
   # there is then a check to see if lattice positions are filled up during the movement process
   # this is sum(stuck)
    fullup=0
    for(pp in (row-1):(row+1))
       {
       for(qq in (col-1):(col+1))
           {
 
           if(pp>=1 && pp<=grid)  # to stop falling off the edge of the grid
              {
              if(qq>=1 && qq<=grid) # ditto
                 {
                 # check available
                 cat(" p ",pp," q ",qq," row ",row," col ",col, "\n")
                 if(sizenum[pp,qq,101]>=100)  # it can only be less than or equal
                   {
                    
                    fullup=fullup+1
                    }
                  }                 
              }
            }
        }
    
    nowheretogo=0
    if(fullup==9)
        {
        nowheretogo=1
        }
    if(nowheretogo==0)
      {      
     
      #  ok there is some space that can be moved to
       nowfullup=0
       while(movers>0)
          {
          trip=1   
          stuck<-matrix(0,nrow=3,ncol=3)   
          while(trip==1 && sum(stuck)<9)  ## two conditions 1 the entity needs to move and sum(stuck) a counter of full up locally
             {
             rr<-round((runif(1)*3)+0.5)-2    
             newr<-row+rr
             cc<-round((runif(1)*3)+0.5)-2
             newc<-col+cc
             stuck[rr,cc]=1        
             if(newr >=1 && newr <=grid)              
               {
               if(newc>=1 && newc <=grid)
                  {
                  if(sizenum[newr,newc,101]<100)
                     {
                     sizenum[newr,newc,101]=sizenum[newr,newc,101]+1
                     place=sizenum[newr,newc,101]
                     sizenum[newr,newc,place]=babycell
                     movers=movers-1
                     trip=0
                    }
                  else
                    {
                    # now full at this point
                    ## we need to allow fill up and then kill
                    stuck[(rr+2),(cc+2)]=1
                    }
                  } #newc
               }# new r
     
             } #while loop
          } # while loop for movers
       }
     else
       {
      # movers float off as there was nowhere to go
      cat("could not find a home for ",movers," \n")
      sizenum[row,col,103]=movers
      movers=0
       }
    }
 }
else
 {
 # just send food back for a bit more diffusion as nothing was able to use it
 sizenum[r,c,102]=avail
 }
return(sizenum)
}



####################################################################
#                                                   PLOT IB OUTPUTS
####################################################################
 plotit <-function(sizenum){        

    # pretty puitvure
    
    squares<-matrix(0,grid,grid)
    i=1
    for(z in 1:grid)
      {
      for(g in 1:grid)
         {
         a=0
         for(kk in 1:100)
           {
           if(sizenum[z,g,kk]>0)
              {
              a=a+1
             
              }
           }
        squares[z,g]=a
 
        i=i+1
    
         }
     }
 

    x <- 10*1:nrow(squares)
    y <- 10*1:ncol(squares)
    filled.contour(x, y, squares, color = terrain.colors)
    #persp(x,y,squares,col=red) 

  }


#persp(x,y,squares,col = terrain.colors) 




##################################################################################################
####################### joint modle runs 
##################################################################################################


##########################################################################################################
                                        DIFFUSION
##########################################################################################################
# this is based on a Soetart function from ReacTran
# there are many more components that could be defined
# key here are the diffusions rates in x,y and z planes
# boundary conditions in 6 faces of box

diffusion3D <- function(t, Y, par) {

  yy    <- array(dim = c(grid, grid, grid), data = Y)  # vector to 3-D array
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
#n=15 # gtid dimensions
dy   <- dx <- dz <- 1   # grid division
Dy   <- Dx <- Dz <- 1   # diffusion coeff, X- and Y-direction

y  <- array(dim = c(grid, grid, grid), data = 0)     # starting concentrations
#
BND   <- matrix(nrow = grid, ncol = grid, 5)      # boundary concentration at upstream point


filled.contour(y[, , 1], main=i,color.palette = terrain.colors)


## parameters and space for IB model############################
sizenum<-array(dim = c(grid, grid, 103), data = 0)  # the array for each grid cell (contains 100 inidivduals)
grate=4         # maximum growth rate per time step
mortality=0.95  # actually survival
maxcell=5       #maximum cell size before division
maintenance=1   # must have this in size to survive a time step


## random populations starting points in grid
pops<-50
for(i in 1:pops)
 {
 rownu<-round((runif(1)*(grid-1))+0.5)
 colnu<-round((runif(1)*(grid-1))+0.5)
 cellno<-round((runif(1)*5)+0.5)
 for(q in 1:cellno)
  {
   sizenum[rownu,colnu,q]=2
  }
 sizenum[rownu,colnu,101]<-cellno
 }


### how long to run
runno=20

popout<-matrix(0,nrow=runno,ncol=4)
popout[,1]=seq(1:runno)

## joint model runs ###########################################

for(i in 1:runno)
  {
  cat("*******************  i ****",i,"\n")
     # STEADY STATE this is the way Dana did it
     diff2.out<-steady.3D(y, func = diffusion3D, parms = NULL, dimens = c(grid, grid, grid),time=0, lrw = 30000000, atol = 1e-10,rtol = 1e-10)
     y <- array(dim=c(grid, grid, grid), data = diff2.out$y)
     # INTEGRATION
     # integration bit more sensitive 
     #diff.out <- ode.3D(y, func = diffusion3D, parms = NULL, dimens = c(grid, grid, grid),times = seq(1,2,0.1), lrw = 30000000, atol = 1e-10,rtol = 1e-10)
     #y <- array(dim = c(grid, grid, grid), data = diff.out[nrow(diff.out), -1])
     
  
    dev.new()
    filled.contour(y[, , 1], main=i,color.palette = terrain.colors)
     # now let the bacteria grow and consume .

    cat("food available at biofilm surface", "\n") 
    # should run this at random through r and c
    randrow<-sample(grid)
    randcol<-sample(grid)
    r=0
    c=0
     for(rr in 1:grid)
        {
        r=randrow[rr]
        for(cc in 1:grid)
           {       
           c=randcol[cc]
            #cat(" well ? ",y[r,c,1]," and whats here ",sizenum[r,c,101],"\n ")
            avail<-y[r,c,1]
            sizenum[r,c,102]<-avail      
 
           #  use the material as individuals ?
            if(sizenum[r,c,101]>0)
               {         
               cat("######food available = ",avail," \n")
               sizenum<-siz_ind(sizenum,grid,r,c, grate,mortality,maxcell,maintenance)
               #food back to diffusion
               cat("food left = ",sizenum[r,c,102]," \n")
                y[r,c,1]=sizenum[r,c,102]
            
               #cat(" ",y[r,c,1],"")
               }
            }
        }
       cat("\n")
       cat("Food after growth ")
       filled.contour(log(y[, , 1]+1), main=i,color.palette = terrain.colors)
       for(pp in 1:grid)
          {
          for(qq in 1:grid)
             {
             cat(" ", y[pp,qq,1], " ")
             }
             cat("\n")
          }
   #get total pop    
   # the cells are growing as a biofilm in each grid of Y at y[1:10,1:10,1] 
  #dev.new()
  popout[i,1]=i
  popout[i,2]=sum(sizenum[,,101])
  popout[i,3]=sum(sizenum[,,102])
  popout[i,4]=sum(sizenum[,,103])
  #plotit(sizenum)

 }


xyplot(popout[,2]+log(popout[,3]+1)+log(popout[,4]+1)~popout[,1], typ="b", ylab="pink nutrients, blue bacteria", xlab="Time")



###########################################################################################

