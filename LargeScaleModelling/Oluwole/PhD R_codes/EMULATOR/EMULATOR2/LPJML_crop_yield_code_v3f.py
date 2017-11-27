class LPJML_crop_yield_code_v3f:

    def __init__(self):
        self._setDurationCalled=False
        self._setConfigCalled=False
        self._aggregated_output=False # default to false, so only override if true

    # optional provision of lat,lon and month (required for run method)
    def setCoordinates(self,lat,lon,month):
        self.lat=lat
        self.lon=lon
        self.month=month
        self._setDurationCalled=True

    # optional control of the configuration via BFG
    def setConfig(self,j,rcp,fertilization,man,coupling,crop_type,path="./",fileName="configuration_option.R"):
        import datetime
        j_options=[0,1,2,3,4,5,6,7,8]
        assert j in j_options, "Error in LPJML_crop_yield_code_v3f:setInputs, value "+str(j)+" is not in "+str(j_options)+"\n"
        rcp_options=["RCP3PD","RCP45","RCP6","RCP85","other"]
        fertilization_options=[0,1]
        man_options=[1,2,3,4,5,6,7]
        coupling_options=["GEMINI","other"]
        crop_type_options=["irrigated","rainfed"]
        f = open(fileName, 'w')
        f.write("# file created by LPJML_crop-yield_code_v3f.py\n")
        f.write("# "+str(datetime.datetime.now())+"\n")
        f.write("# Any changes made to this file will be overwritten when running the python script\n")
        self._add_option(f,"path",path)
        self._add_option(f,"j",j,options=j_options)
        self._add_option(f,"rcp",rcp,options=rcp_options)
        self._add_option(f,"fertilization",fertilization,options=fertilization_options)
        self._add_option(f,"man",man,options=man_options)
        self._add_option(f,"coupling",coupling,options=coupling_options)
        if coupling=="GEMINI" : self._aggregated_output=True
        self._add_option(f,"type",crop_type,options=crop_type_options)
        f.close()
        self._setConfigCalled=True

    # input data is passed via files and we synchronise using BFG
    def runFileSync(self,go):
        self.go=self._setValue(go)
        assert self.go==True,"Error in LPJML_crop_yield_code_v3f:runFileSync() expecting the 'go' argument to be set to True"
        self._runLPJML()
        return self._getResults()

    # input data is passed via BFG
    def run(self,cloud,precipitation,temperature,wet_days):
        assert self._setDurationCalled, "Error in LPJML_crop_yield_code_v3f:run(). _setDuration() must be called before this function"
        self._setInputs(cloud,precipitation,temperature,wet_days)
        self._runLPJML()
        return self._getResults()

    # private methods

    def _setInputs(self,cloud,precipitation,temperature,wet_days):
        self._netcdf_write(cloud,"cld.nc","CLD_MONTHLY","percent")
        self._netcdf_write(precipitation,"pre.nc","PRECIP_MONTHLY","mm")
        self._netcdf_write(temperature,"tmp.nc","TMP_MONTHLY","celsius degrees")
        self._netcdf_write(wet_days,"wet.nc","WET_MONTHLY","days")

    def _netcdf_write(self,data,fileName,dataName,dataUnits):

        #from Scientific.IO.NetCDF import NetCDFFile
        from scipy.io import netcdf

        # note: we overwrite any existing file
        assert len(data.shape)==3,"Error, cloud expected to be 3D"
        nlon=data.shape[0]
        nlat=data.shape[1]
        nmonth=data.shape[2]
        print "python:netcdf_write:data "+dataName+" nlat="+str(nlat)+" nlon="+str(nlon)+" nmonth="+str(nmonth)

        file = netcdf.netcdf_file(fileName, 'w')
        #file = NetCDFFile(fileName, 'w')

        file.createDimension('LONGITUDE',nlon)
        file.createDimension('LATITUDE',nlat)
        file.createDimension('MONTH',nmonth)

        assert len(self.lat)==nlat,"Error, lat array is not the same size as nlat"
        dims = 'LATITUDE',
        dataVar = file.createVariable('LATITUDE', 'f', dims)
        #dataVar.assignValue(self.lat)
        dataVar[:]=self.lat
        dataVar.units='degrees'

        assert len(self.lon)==nlon,"Error, lon array is not the same size as nlon"
        dims = 'LONGITUDE',
        dataVar = file.createVariable('LONGITUDE', 'f', dims)
        #dataVar.assignValue(self.lon)
        dataVar[:]=self.lon
        dataVar.units='degrees'

        assert len(self.month)==nmonth,"Error, month array is not the same size as nmonth"
        dims = 'MONTH',
        dataVar = file.createVariable('MONTH', 'f', dims)
        #dataVar.assignValue(self.month)
        dataVar[:]=self.month
        dataVar.units='years'

        dims = ('LONGITUDE','LATITUDE','MONTH')
        dataVar = file.createVariable(dataName, 'f', dims)
        #dataVar.assignValue(data)
        dataVar[:]=data
        dataVar.units=dataUnits

        file.close()


    def _runLPJML(self,skip=False):
        if skip : return
        import subprocess
        import os
        import sys
        # Call the R script with our arguments
        infile = open("LPJML_crop-yield_code_v3f.R", 'r')
        outfile = open("out.txt", 'w')
        args = ["R", "--no-restore", "--no-save"]
        proc = subprocess.Popen(args, stdin=infile, stdout=outfile, env=os.environ, close_fds=True)
        proc.wait()

        # check for success
        if proc.returncode:
            print("[BFG] Subprocess \"{0}\" in module {1} returned error {2}, aborting.".format(" ".join(args), __name__, proc.returncode))
            sys.exit(proc.returncode)

    def _add_option(self,fileHandle,name,value,options=None):
        if options is not None:
            assert value in options, "Error in LPJML_crop_yield_code_v3f:add_option, value "+str(value)+" is not in "+str(options)+"\n"
            fileHandle.write("# "+name+" options "+str(options)+"\n")
        if isinstance(value, str): # python 3 assumption
            fileHandle.write(name+"=\""+value+"\"\n")
        else:
            fileHandle.write(name+"="+str(value)+"\n")

    def _isAggregatedOutput(self,fileName="configuration_option.R"):
        # read config file to find out if coupling is set to GEMINI or not
        import re
        f = open(fileName, 'r')
        p = re.compile('^[ \t]*coupling')
        result=None
        for i, line in enumerate(f):
            if p.match(line):
              result=line.split("\"")[1]
        f.close()
	assert result is not None, "Error in _isAggregatedOutput, 'coupling' option not found in file "+fileName
	if result=="GEMINI":
            return True
	elif result=="other":
            return False
        else:
	    assert False, "Error in _isAggregatedOutput, expected 'coupling' option in file "+ fileName+" to be one of [GEMINI,other] but found '"+result+"'"

    def _getResults(self):
        lat_csv=[];lon_csv=[];cereal_csv=[];rice_csv=[];maize_csv=[];oil_csv=[]
        if not self._setConfigCalled : self._aggregated_output=self._isAggregatedOutput()
        if self._aggregated_output:
	    return self._csv_read()
        else: # netcdf has been written
	    # assume nothing needs to be returned for the timebeing.
            # We will need to return both the aggregated and standard output if standard
            # output is needed for coupling.
            #lat_ncf,lon_ncf,cereal_ncf,rice_ncf,maize_ncf,oil_ncf=self._netcdf_read()
            return [],[],[],[],[],[]

    def _netcdf_read(self):
        assert False, "_getncfResults() not yet implemented"

    def _csv_read(self):
        import os
        import csv
        # extract the required results and return
        fileSearch=os.path.join('output_result','*.csv')
        csvFileNames=glob.glob(fileSearch)
        assert len(csvFileNames)==1, "Error in LPJML_crop_yield_code_v3f:_getResults(), expecting one output file to match "+fileSearch+" but found "+str(len(csvFileNames))+" : "+str(csvFileNames)
        print "python:_getResults:csv read: "+str(csvFileNames)
        fileName=csvFileNames[0]
        reader = csv.reader(open(fileName,"rb"),delimiter=",")
        #format expected is : lat, lon, cereal, rice, maize, oil
        lat=[];lon=[];cereal=[];rice=[];maize=[];oil=[]
        first=True
        for row in reader:
            if first:
                # skip the first row as it contains the titles
                first=False
            else:
                lat.append(row[0])
                lon.append(row[1])
                cereal.append(row[2])
                rice.append(row[3])
                maize.append(row[4])
                oil.append(row[5])
        return lat,lon,cereal,rice,maize,oil

    # hack as python booleans are being returned as ints by the python layer in BFG
    def _setValue(self,origValue):
        if type(origValue) is bool:
            newValue=origValue
        elif type(origValue) is int:
            if origValue==0:
                newValue=False
            else:
                newValue=True
        else:
            raise Exception, "Error"
        return newValue

if __name__ == "__main__":

    # create model
    myLPJML_crop_em = LPJML_crop_yield_code_v3f()

    # run model Sudipta style - assumes input and config files pre-exist
    lat,lon,cereal,rice,maize,oil=myLPJML_crop_em.runFileSync(True)

    # Do something with the output
    print "len(lat)="+str(len(lat))
    print "len(lon)="+str(len(lon))
    print "len(cereal)="+str(len(cereal))
    print "len(rice)="+str(len(rice))
    print "len(maize)="+str(len(maize))
    print "len(oil)="+str(len(oil))

    # Can not continue as we need to remove the pre-existing netcdf input files.
    # The setConfig() method also overwrites any pre-existing config file.
    exit(0)

    # run model with input data
    import numpy
    # set up input data
    nlat=280
    nlon=720
    nmonth=1200
    lat=numpy.linspace(-55.75,83.75,nlat) # how do I make these float32?
    lon=numpy.linspace(-179.75,179.75,nlon)
    month=numpy.linspace(2001.0,2100.917,nmonth)
    cld=numpy.zeros((nlon,nlat,nmonth),numpy.float32)
    pre=numpy.zeros((nlon,nlat,nmonth),numpy.float32)
    tmp=numpy.zeros((nlon,nlat,nmonth),numpy.float32)
    wet=numpy.zeros((nlon,nlat,nmonth),numpy.float32)
    # initialise input coordinates - required before run()
    myLPJML_crop_em.setCoordinates(lat,lon,month)
    # run model standard BFG style - creates input files
    #lat,lon,cereal,rice,maize,oil=myLPJML_crop_em.run(cld,pre,tmp,wet)

    # Do something with the output
    #print "len(lat)="+str(len(lat))
    #print "len(lon)="+str(len(lon))
    #print "len(cereal)="+str(len(cereal))
    #print "len(rice)="+str(len(rice))
    #print "len(maize)="+str(len(maize))
    #print "len(oil)="+str(len(oil))

    # set configuration values BFG style - overwrites config file
    j=2;rcp="RCP6";fertilization=1;man=5;coupling="other";crop_type="rainfed"
    myLPJML_crop_em.setConfig(j,rcp,fertilization,man,coupling,crop_type)
    # run model standard BFG style - you could run Sudipta style here instead i.e. setConfig() can be used with runFileSync()
    lat,lon,cereal,rice,maize,oil=myLPJML_crop_em.run(cld,pre,tmp,wet)

    # Do something with the output
    print "len(lat)="+str(len(lat))
    print "len(lon)="+str(len(lon))
    print "len(cereal)="+str(len(cereal))
    print "len(rice)="+str(len(rice))
    print "len(maize)="+str(len(maize))
    print "len(oil)="+str(len(oil))

    print "complete"
