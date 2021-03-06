#!/usr/bin/python
import argparse
import sys
import codecs
import numpy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy.optimize as opt
from math import log10, floor 
import ntpath
import os


#Parameters to set:
TgRANGEhigh=80 #top temperature to examine. Default is 80
TgRANGElow=20#lowest temperature to examine. Default is 20
degforlin=2 #window in degrees up or down from midpoint to calculate linear regions. Default is 2
degspace=5 #space in degrees between linear regions. Default is 5
smooth=20 # value for smoothing function, 1 is no smoothing. Higher averages nearby points. Default is 20.
deviation=0.1 # devation from baseline (preline) towards end line (postline), 0.1 = 10% shift towards post line. Default is 0.1
#tempshift=1 #used to scale temperature to better solve

#round to 3 sig figs
def round_sig(x, sig=3):
    return round(x, sig-int(floor(log10(abs(x))))-1)

def path_leaf(path):
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)

#function to test the fit of lines
def fit(params,midtemp,temprange,Temp=None,hflow=None,fixedm=None,): #((m,b),midtemp,temprange), Returns total difference
    b = params[1] #second parameter is "b"
    if fixedm:
        m=fixedm
    else:
        m = params[0] # first parameter is m, unless given as fixedm
    comparerange=[i for i,x in enumerate(Temp) if midtemp-temprange<x<midtemp+temprange] # list with the positions in the temperature window for the low temp line
    line = m*numpy.array(Temp)+b # calculate line (y=mx+b) store y values as line
    totaldiff= 0
    for i in comparerange:
        totaldiff+=(line[i]-hflow[i])**2
    totaldiff=100*totaldiff/len(comparerange)
    return totaldiff

#find linear regions
def linehunt(midstart,mbbounds,searchspace=degspace,temprange=degforlin,Temp=None,hflow=None,m=None):
    fits=[]
    for i in range(0,(1+searchspace)): #iterate in 1 degree steps over the searchspace (+1 to include max value)
        midtemptest=midstart+i
        parameters=(0,0) # best guesses, could add at argument
        if m:
            attempt = opt.minimize(fit, parameters, args=(midtemptest,temprange,Temp,hflow,m), method='SLSQP', bounds=mbbounds, options={'maxiter': 1e3})
            #if attempt.x[0]!=m:
            #    print "solved m is: "+str(attempt.x[0])+", but fixed m is: "+str(m)
        else:
            attempt = opt.minimize(fit, parameters, args=(midtemptest,temprange,Temp,hflow), method='SLSQP', bounds=mbbounds, options={'maxiter': 1e3})
        fits.append(attempt.fun) # make list of best fit scoreing
        parameters=(attempt.x[0],attempt.x[1])#update parameters to allow next search to be faster
    bestfit=min(fits)
    besti=[i for i, x in enumerate(fits) if x==bestfit]
    midtemp=midstart+besti

    #do a better, more honed in fitting once midtemp is identified
    if m:
        result = opt.minimize(fit, parameters, args=(midtemp,temprange,Temp,hflow,m), method='SLSQP', bounds=mbbounds, options={'maxiter': 1e6})
    else:
        result = opt.minimize(fit, parameters, args=(midtemp,temprange,Temp,hflow), method='SLSQP', bounds=mbbounds, options={'maxiter': 1e6})
    print "best fit: "+str(round_sig(result.fun,2)),
    print "fit range "+str(round_sig((midtemp-temprange)[0]))+" to "+str(round_sig((midtemp+temprange)[0]))+" deg C. Equation:",
    if m:
        print "y = "+str(round_sig(x))+"x  +"+str(round_sig(result.x[1]))
        return(midtemp,m,result.x[1],temprange) # returns midtemp of line, "m" and, solved "b", and temp range (up and down from midtemp)
    else:    
        print "y = "+str(round_sig(result.x[0]))+"x + "+str(round_sig(result.x[1]))
        return(midtemp,result.x[0],result.x[1],temprange)# returns midtemp of line, solved "m" and, solved "b", and temp range (up and down from midtemp)

#main function to analyze a DSC file
def Analyze(filename,graphs=True):
    basename=path_leaf(filename)
    print "opening file \""+basename+"\"",
    data=[]
    lines=[]
    
    try:
        f = codecs.open(filename, 'r', 'utf-8') # open the file given, codecs reads it as utf-8
        temp = f.read()#try reading to check encodings
    except UnicodeError:
        print "file not standard unicode, trying to read as utf-16"
        try:
            f = codecs.open(filename, 'r', 'utf-16') # open the file given, codecs reads it as utf-16
            contents = f.read()
        except UnicodeError:
            print "file not utf-18 either"
            sys.exit(1)
        else:
            f = codecs.open(filename, 'r', 'utf-16') # Error unless reopened file in an else statement. It's like it closes it after the try
    else:
        f = codecs.open(filename, 'r', 'utf-8')  # Error unless reopened file in an else statement. It's like it closes it after the try
    
    #get rid of whitespace & newline charcters
    for line in f:
        #print line
        lines.append(line.strip('\n')) #this only works if encoding is correct
    #Find the Start of data "StartOfData"
    keep=0
    key= "StartOfData" #text that indicates start of data
    for line in lines:
        if keep: #after start of data is found, keep is 1 and lines are added to "data"
            data.append(line)
        if key in line: #if the key is found, switch keep to 1 to begin saving lines
            print "-->> found data"
            keep=1
        
    raw = numpy.genfromtxt(data, delimiter = "\t") #convert list of tab deliminited "data" into numby array
    '''
    Positions in raw array:
    Sig1	Time (min)
    Sig2	Temperature (C)
    Sig3	Heat Flow (mW)
    Sig4	Heat Capacity (mJ/C)
    Sig5	Sample Purge Flow (mL/min)
    Sig6	LNCS Pressure (kPa gauge)
    '''
    Temp=raw[:,1]
    hflow=raw[:,2]
    
    #calculate derivative set
    hflowd=hflow.copy() #make list of derivative values by coping hflow list
    
    for i in range(len(hflowd[:])):
        if i < smooth: # set initial points that can't be smoothed to 0
            hflowd[i]=0
            #print 'start'
        elif i> len(hflowd[smooth-1:]): # set last points that can't be smoothed to 0
            hflowd[i]=0
            #print "end"
        else: #calculate smoothed derivative for all other points
            dif1=sum(hflow[i-smooth:i]) #sum "smooth" number of points before
            tempdif1=sum(Temp[i-smooth:i]) # sum temperatures for each point
            dif2=sum(hflow[i:i+smooth]) # sum "smooth" number of points after
            tempdif2=sum(Temp[i:i+smooth]) # sum temperatures for each point
            hflowd[i]=(dif2-dif1)*(tempdif2-tempdif1)**-1 #set values to smoothed derivative difference calculated between [smooth] points that are [smooth] apart divided by sum of temps normalize.
    
    Tg=0 #max dir
    window=[i for i,x in enumerate(Temp) if TgRANGElow<x<TgRANGEhigh] # list with the positions in the temperature window
    Tgflow=max(hflowd[window[:]]) # Tgflow set to maximum change in heat flow
    Tgpos= [i for i,x in enumerate(hflowd) if x==Tgflow] # Tgpos set to position in list for max change in heat flow
    Tg= Temp[Tgpos[0]] # Tg set to temperature at maximum change in heat flow, i.e. steepest point in curve
    print "Tg (maximum slope): "+str(round(Tg,1))+" deg C"
    #define regions to determine linear reagions
    transrange=[i for i,x in enumerate(Temp) if Tg-(0.5*degforlin)<x<Tg+(0.5*degforlin)] # list with the positions in the temperature window for the low temp line
    
    
    #point to start looking for linear ranges
    premid=Tg-2*degforlin-degspace#start looking
    postmid=Tg+2*degforlin
    
    #fit lines
    bnd = (-10, 10)
    prebnds=(bnd,bnd)#,(premid-degspace,premid),(degforlin,degforlin))
    transbnds=(bnd,bnd)#,(Tg,Tg),(degforlin,degforlin))
    postbnds=(bnd,bnd)#,(postmid,postmid+degspace),(degforlin,degforlin))
    
    
    print "FITTING GLASS LINE"
    preline=linehunt(premid,prebnds,Temp=Temp,hflow=hflow)
    premid, prem, preb, pretemprange = preline
    
    print "FITTING LIQUID LINE"
    postline=linehunt(postmid,postbnds,Temp=Temp,hflow=hflow)
    postmid, postm, postb, posttemprange = postline
    
    print "FITTING TRANSITION LINE"
    transline=linehunt(Tg,transbnds,searchspace=0,Temp=Temp,hflow=hflow,m=Tgflow)
    transmid, transm, transb, transtemprange = transline
    #sanity check
    if transm != Tgflow:
        print "we have a problem. transm is: "+str(transm)+" while Tgflow is: "+str(Tgflow)
    
    #ranges used for pre and post line fitting
    prerange=[i for i,x in enumerate(Temp) if premid-pretemprange<x<premid+pretemprange] # list with the positions in the temperature window for the low temp line
    postrange=[i for i,x in enumerate(Temp) if postmid-posttemprange<x<postmid+posttemprange] # list with the positions in the temperature window for the low temp line
    
    #generate datasets from best fit lines for pre, trans, and post ranges
    pre = prem*numpy.array(Temp)+preb # calculate line pre (y=mx+b) store y values as pre
    post = postm*numpy.array(Temp)+postb # calculate line pre (y=mx+b) store y values as pre
    trans = transm*numpy.array(Temp)+transb # calculate line pre (y=mx+b) store y values as pre
    
    deviationstart=[i for i,x in enumerate(Temp) if premid<x<Tg] # list with the positions in the temperature window to look for devations from pre line
    print "" 
    print"****  FILE: "+basename+"  *****"
    print"****  Tg is: "+str(round(Tg,1))+" degrees, ",
    Onset=None
    End=0 
    #find the onset as defined by 10% of the way between pre line and post line
    for i in deviationstart:
        if not Onset:
            if hflow[i]>((deviation*(postb+postm*Temp[i]))+((1-deviation)*(preb+prem*Temp[i]))):
                Onset=Temp[i]
    if Onset: 
        print "Glass Transition Onset: "+str(round(Onset,1))+" deg C  ****"
    else: print" Onset not found  ****"
    print ""
    
    if graphs:
        ##### Draw Figure
        plt.figure() # derivative plot
        plt.plot(Temp, hflowd, label = "derivative")
        plt.axvline(x=Tg,color='burlywood',label="Tg")
        plt.axvline(x=Onset,color="orange",label="Tg_onset")
        #plt.axvspan(Onset, End, facecolor='#2ca02c', alpha=0.2)
        #plt.axvline(x=Onset, color='chartreuse',label="Onset")
        plt.legend(loc=2)
        
        # x-axis label
        plt.xlabel('Temp deg C')
        # frequency label
        plt.ylabel('heat flow derivative (mW/C)')
        # plot title
        plt.title(basename+" derivative")
        
        plt.xlim(TgRANGElow,TgRANGEhigh)
        plt.ylim((min(hflowd[window[:]])-0.01),(max(hflowd[window[:]])+0.01)) # y window set to maximum range +/- 0.01
        plt.savefig(basename+"_derivative.png")
        plt.close()
        
        plt.figure() # main plot
        #remove unwanted portions of lines
        pre[:prerange[0]]=numpy.nan #remove pre line up through the begining of pre range
        pre[Tgpos[0]:]=numpy.nan #remove pre line past Tg
        trans[:0.5*(transrange[1]+prerange[-1])]=numpy.nan #remove trans line up through half-way between the end of pre range and the start of trans range
        trans[0.5*(postrange[0]+transrange[-1]):]=numpy.nan #remove trans line past half-way between the end of trans range and the start of post range
        post[:Tgpos[0]]=numpy.nan # remove post line up through Tg
        post[postrange[-1]:] =numpy.nan# remove post line after post range
        matplotlib.pyplot.plot(Temp,hflow,label = "data")
        plt.plot(Temp, pre, label = "glass",linestyle='dashed', lw=2)
        plt.plot(Temp, post, label = "liquid", linestyle='dashed', lw=2)
        plt.plot(Temp, trans, label = "transition", linestyle='dashed', lw=2)
        
        plt.axvline(x=Tg,color='burlywood',label="Tg")
        plt.axvline(x=Onset, color="orange", label="Tg_onset")
        plt.axvspan(Temp[prerange[0]], Temp[prerange[-1]], facecolor='#0c002c', alpha=0.2)
        plt.axvspan(Temp[postrange[0]], Temp[postrange[-1]], facecolor='#0ca000', alpha=0.2)
        plt.legend(loc=2)
        
        # x-axis label
        plt.xlabel('Temp deg C')
        # frequency label
        plt.ylabel('heat flow (mW)')
        # plot title
        plt.title(basename)
        
        #plt.xlim((-10, 100))
        plt.xlim(TgRANGElow,TgRANGEhigh)
        plt.ylim((min(hflow[window[:]])-0.1),(max(hflow[window[:]])+0.1))  # y window set to maximum range +/- 0.1
        plt.savefig(basename+"_plot.png")
        plt.close()

#READING DATA IN
parser = argparse.ArgumentParser(description='File',prog='Analyze.py') #read the flags
parser.add_argument('-f', '--file', metavar='FILENAME.txt',dest='FILENAME',help='txt data file to be analyzed') # read flag for file into "FILENAME" variable
args = parser.parse_args()

if args.FILENAME==None:
    print "no file specified, analyzing files in ToProcess folder"
    os.listdir("./ToProcess")
    filelist = [f for f in os.listdir("./ToProcess") if os.path.isfile(os.path.join("./ToProcess", f))]
    plotnumber=0
    for file in filelist:
        Analyze("./ToProcess/"+file)
    sys.exit()
elif not os.path.isfile(args.FILENAME):
    sys.exit(args.FILENAME+" does not exist. Please specify a DSC file, or leave blank and place file(s) to be analyzed in /ToProcess")
else:
    Analyze(args.FILENAME)
    sys.exit()