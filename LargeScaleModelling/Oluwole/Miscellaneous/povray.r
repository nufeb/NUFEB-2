## R-POVRAY
require(grDevices)
# Scene class
Scene=function(background=Colour(1,1,1)){
  self=list(
    contents=c(paste("background{",background$format(light=T),"}\n",sep="")),
 #gaja="C:\\Users\\olu\\Desktop\\NEW_LAMMPS_RESULTS",
gaja=getwd(),
  #filepath=tempfile(pattern="plt",tmpdir=gaja,fileext=".pov"),
filepath=paste(paste(gaja,"plot",sep="\\"), formatC(i, digits = 1, flag = "0"), ".pov", sep = ""),
    add=function(text){
      self$contents=c(self$contents,text)
    },
    save=function(path=self$filepath){
      self$filepath = path
      f=file(path,"w")
      for(line in self$contents){
        cat(line,file=f,sep="\n")
      }
       close(f)
    },
    push=function(object){
      self$contents=c(self$contents,object$format())
    },
    render=function(w=640,h=400,aa=F){
      self$save()
      if(aa){
        aa=" +A0.1 "
      }else{
        aa=" "
      }
      if(Sys.info()['sysname'] %in% c("Darwin","Linux") ){
        POVLOC = Sys.which("povray")
        command=paste(" +I",self$filepath,aa," +W",w," +H",h," +P -V",sep="")
      }else{
        POVLOC="C:\\Program Files\\POV-Ray\\v3.7 RC7\\bin\\pvengine.exe"
command2= paste(paste(gaja,"plot",sep="\\"), formatC(i, digits = 1, flag = "0"), ".pov", sep = "")
        command=paste("",command2,aa," +W",w," +H",h," +P -V /EXIT",sep="")
      }
#command2 <- paste("plot", formatC(i, digits = 1, flag = "0"), "", sep = "")
      #out=system2(POVLOC,args=command)
out=system2(POVLOC,args=command)
      out
    }
    )
  self <- list2env(self)
  class(self) <- "Scene"
  return(self) 
}  

# Light class
Light=function(location=c(x,y,z),col=NA){
  self=list(
    location=location,
    col=col,
    format=function(){
      x=self$location[1]
      y=self$location[2]
      z=self$location[3]
      text=paste("light_source {  <",x,",",y,",",z,">",sep="")
      if(is.environment(col)){
        text=paste(text,col$format(light=T))
      }
      text=paste(text,"}")
      text
    }
  )
  self <- list2env(self)
  class(self) <- "Light"
  return(self) 
}

# Camera class
Camera=function(location=c(x,y,z),looking=c(x2,y2,z2)){
  self=list(
    location=location,
    looking=looking,
    format=function(){
      x=self$location[1]
      y=self$location[2]
      z=self$location[3]
      x2=self$looking[1]
      y2=self$looking[2]
      z2=self$looking[3]
      text=paste("camera { location <",x,",",y,",",z,"> look_at <",x2,",",y2,",",z2,"> right x*image_width/image_height}",sep="")
      text
    }
  )
  self <- list2env(self)
  class(self) <- "Camera"
  return(self) 
}

Colour=Color=function(r=1,g=1,b=1,t=NA,f=NA){
  self=list(
    r=r,
    g=g,
    b=b,
    t=t,
    f=f,
    format=function(light=F){
      output=""
      if(light==F){
        output="texture{pigment{"
      }
      if((is.na(self$t) & is.na(self$f)) | (!is.na(self$t) & !is.na(self$f))){
        output = paste(output,"color rgb<",self$r,",",self$g,",",self$b,">",sep="")
      }
      if(!is.na(t) & is.na(f)){
        output = paste(output,"color rgbt<",self$r,",",self$g,",",self$b,",",self$t,">",sep="")
      }
      if(is.na(t) & !is.na(f)){
        output = paste(output,"color rgbf<",self$r,",",self$g,",",self$b,",",self$f,">",sep="")
      }
      if(light==F){
        output=paste(output,"}}")
      }
      output
    }
    )
  self <- list2env(self)
  class(self) <- "Colour"
  return(self) 
}

Texture=function(tex){
  self=list(
    tex=tex,
    format=function(){
      self$tex
    }
  )
  self <- list2env(self)
  class(self) <- "Texture"
  return(self) 
}

Interior=function(interior){
  self=list(
    interior=interior,
    format=function(){
      paste("interior{",self$interior,"}")
    }
  )
  self <- list2env(self)
  class(self) <- "Interior"
  return(self) 
}

gen_texture=function(col,tex,interior){
  text=""
  if(is.environment(col) & !is.environment(tex)){
    text=paste(text,col$format(light=F))
  }
  if(is.environment(tex) & !is.environment(col)){
    text=paste(text,"texture{",tex$format(),"}")
  }
  if(is.environment(tex) & is.environment(col)){
    text=paste(text,"texture{",tex$format()," pigment{",col$format(light=T),"}}")
  }
  if(is.environment(interior)){
    text=paste(text,interior$format())
  }
  text
}

Cylinder=function(start=c(x1,y1,z1),end=c(x2,y2,z2),radius=1,col=NA,tex=NA,interior=NA){
  self=list(
    start=start,
    end=end,
    col=col,
    tex=tex,
    interior=interior,
    radius=radius,
    transform=c(),
    translate=function(translate=c(x,y,z)){
      self$transform=c(self$transform,paste("translate{<",translate[1],",",translate[2],",",translate[3],">}"))
    },
    rotate=function(rotate=c(x,y,z)){
      self$transform=c(self$transform,paste("rotate{<",rotate[1],",",rotate[2],",",rotate[3],">}"))
    },    
    scale=function(scale=c(x,y,z)){
      if(length(scale)==3){
        self$transform=c(self$transform,paste("scale{<",scale[1],",",scale[2],",",scale[3],">}"))
      }else{
        self$transform=c(self$transform,paste("scale{",scale[1],"}"))
      }    
    },
    format=function(){
      x1=self$start[1]
      y1=self$start[2]
      z1=self$start[3]
      x2=self$end[1]
      y2=self$end[2]
      z2=self$end[3]
      text=""
      text=paste(text,"cylinder{",sep="")
      text=paste(text,"<",x1,",",y1,",",z1,">,<",x2,",",y2,",",z2,">,",self$radius,sep="")
      text=paste(text,gen_texture(self$col,self$tex,self$interior),sep=" ")
      tl=paste(self$transform,sep=" ",collapse=" ")
      text=paste(text,tl,sep=" ")
      text=paste(text,"} ",sep=" ")
      text
    }
    
    )
  self <- list2env(self)
  class(self) <- "Cylinder"
  return(self) 
}

Cone=function(start=c(x1,y1,z1),end=c(x2,y2,z2),radius=1,radius2=1,col=NA,tex=NA,interior=NA){
  self=list(
    start=start,
    end=end,
    col=col,
    tex=tex,
    interior=interior,
    radius=radius,
    radius2=radius2,
    transform=c(),
    translate=function(translate=c(x,y,z)){
      self$transform=c(self$transform,paste("translate{<",translate[1],",",translate[2],",",translate[3],">}"))
    },
    rotate=function(rotate=c(x,y,z)){
      self$transform=c(self$transform,paste("rotate{<",rotate[1],",",rotate[2],",",rotate[3],">}"))
    },    
    scale=function(scale=c(x,y,z)){
      if(length(scale)==3){
        self$transform=c(self$transform,paste("scale{<",scale[1],",",scale[2],",",scale[3],">}"))
      }else{
        self$transform=c(self$transform,paste("scale{",scale[1],"}"))
      }    
    },
    format=function(){
      x1=self$start[1]
      y1=self$start[2]
      z1=self$start[3]
      x2=self$end[1]
      y2=self$end[2]
      z2=self$end[3]
      text=""
      text=paste(text,"cone{",sep="")
      text=paste(text,"<",x1,",",y1,",",z1,">,",self$radius,",<",x2,",",y2,",",z2,">,",self$radius2,sep="")
      text=paste(text,gen_texture(self$col,self$tex,self$interior),sep=" ")
      tl=paste(self$transform,sep=" ",collapse=" ")
      text=paste(text,tl,sep=" ")
      text=paste(text,"}",sep=" ")
      text
    }
    
  )
  self <- list2env(self)
  class(self) <- "Cylinder"
  return(self) 
}

Box=function(start=c(x1,y1,z1),end=c(x2,y2,z2),col=NA,tex=NA,interior=NA){
  self=list(
    start=start,
    end=end,
    col=col,
    tex=tex,
    interior=interior,
    transform=c(),
    translate=function(translate=c(x,y,z)){
      self$transform=c(self$transform,paste("translate{<",translate[1],",",translate[2],",",translate[3],">}"))
    },
    rotate=function(rotate=c(x,y,z)){
      self$transform=c(self$transform,paste("rotate{<",rotate[1],",",rotate[2],",",rotate[3],">}"))
    },    
    scale=function(scale=c(x,y,z)){
      if(length(scale)==3){
        self$transform=c(self$transform,paste("scale{<",scale[1],",",scale[2],",",scale[3],">}"))
      }else{
        self$transform=c(self$transform,paste("scale{",scale[1],"}"))
      }    
    },
    format=function(){
      x1=self$start[1]
      y1=self$start[2]
      z1=self$start[3]
      x2=self$end[1]
      y2=self$end[2]
      z2=self$end[3]
      text=""
      text=paste(text,"box{",sep="")
      text=paste(text,"<",x1,",",y1,",",z1,">,<",x2,",",y2,",",z2,"> ",sep="")
      text=paste(text,gen_texture(self$col,self$tex,self$interior),sep=" ")
      tl=paste(self$transform,sep=" ",collapse=" ")
      text=paste(text,tl,sep=" ")
      text=paste(text,"}",sep=" ")
      text
    }
    
  )
  self <- list2env(self)
  class(self) <- "Box"
  return(self) 
}

Sphere=function(centre=c(x1,y1,z1),radius=1,col=NA,tex=NA,interior=NA){
  self=list(
    centre=centre,
    col=col,
    tex=tex,
    radius=radius,
    interior=interior,
    transform=c(),
    translate=function(translate=c(x,y,z)){
      self$transform=c(self$transform,paste("translate{<",translate[1],",",translate[2],",",translate[3],">}"))
    },
    rotate=function(rotate=c(x,y,z)){
      self$transform=c(self$transform,paste("rotate{<",rotate[1],",",rotate[2],",",rotate[3],">}"))
    },    
    scale=function(scale=c(x,y,z)){
      if(length(scale)==3){
        self$transform=c(self$transform,paste("scale{<",scale[1],",",scale[2],",",scale[3],">}"))
      }else{
        self$transform=c(self$transform,paste("scale{",scale[1],"}"))
      }    
    },
    format=function(){
      x1=self$centre[1]
      y1=self$centre[2]
      z1=self$centre[3]
      text=""
      text=paste(text,"sphere{",sep="")
      text=paste(text,"<",x1,",",y1,",",z1,">,",self$radius,sep="")
      text=paste(text,gen_texture(self$col,self$tex,self$interior),sep=" ")
      tl=paste(self$transform,sep=" ",collapse=" ")
      text=paste(text,tl,sep=" ")
      text=paste(text,"}",sep=" ")
      text
    }
    
  )
  self <- list2env(self)
  class(self) <- "Sphere"
  return(self) 
}

norm = function(x){
  mn = min(x)
  mx = max(x)
  x=x-mn
  x=x/(mx-mn)
  x
}

povplot=function(x,y,z,obj="Sphere",col=Colour(c(1,0,0)),f=NA,t=NA,tex=NA,size=0.1){
  require(grDevices)
  if(is.environment(col)){
    envcol=T
  }else{
    envcol=F
    if(length(col) == length(x)){
      varcol=T
    }else{
      varcol=F
    }
    colmat=grDevices::col2rgb(col)
  }
  outlist=list()
  x=norm(x)
  y=norm(y)
  z=norm(z)
  for(idx in seq_along(x)){
    if(envcol==T){
      thiscol=col
    }else{
      if(varcol==T){
        place=idx
      }else{
        place=1
      }
      r=as.numeric(colmat[1,place])
      g=as.numeric(colmat[2,place])
      b=as.numeric(colmat[3,place])
      print(c(r,g,b))
      thiscol=Colour(r,g,b,f,t)
    }
    if(obj == "Sphere"){
      print(thiscol$r)
      print(tex)
      os = Sphere(c(x[idx],y[idx],z[idx]),size,col=thiscol,tex=tex)
      outlist=c(outlist,os)
    }
  }
  outlist
}

source("textures.r")