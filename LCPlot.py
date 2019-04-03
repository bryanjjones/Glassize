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
import csv

#Parameters to set:
ElutionRANGEhigh=75 #end of elution to examine. 
ElutionRANGElow=20#beginning of elution to examine. 
Adducts_start=50
Adducts_end=58
IgG_start=58
IgG_end=72
adducts_color='#906000'
IgG_color='#006010'
colors=('000000','000000','000000','000000','000000','000000','000000','000000','000000','000000','000000','000000')
#colors=('#ff0000','#600000','#606000','#60ff00')
#round to 3 sig figs
def round_sig(x, sig=3):
    return round(x, sig-int(floor(log10(abs(x))))-1)

def path_leaf(path):
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)

#test function
def Test(filename,graphs=True):
    basename=path_leaf(filename)
    print "opening file \""+basename+"\""
    data=[]
    lines=[]
    mAU=[]
    mL=[]
    with open(filename, 'rb') as csvfile:
        data = list(csv.reader(csvfile))
    name=(data[0][0])
    for i in range(len(data[3:])):
        if data[i+3][0]:
            mL.append(float(data[i+3][0]))
            mAU.append(float(data[i+3][1]))
    for i in range(len(mL[:])):
        try:
            x=float(mL[i])
        except:
            print "cant turn this into a float"
            print mL[i]
            print i
    return(mL,mAU,name)
    #sys.exit()
def Plot(X,Y,seriesname,color):   
    #line_color=''
    adductsX=[]
    adductsY=[]
    IgGX=[]
    IgGY=[]
    p0=matplotlib.pyplot.plot(X,Y,label = None, color=color, antialiased=True)
    for i in range(len(X[:])):
        if Adducts_start <= X[i] <= Adducts_end:
            adductsX.append(X[i+1] )
            adductsY.append(Y[i+1])
        if IgG_start < X[i] <= IgG_end:
            IgGX.append(X[i])
            IgGY.append(Y[i])
    plt.fill_between(adductsX, adductsY, alpha=.3, color=adducts_color, antialiased=True)
    plt.fill_between(IgGX, IgGY, alpha=.3, color=IgG_color, antialiased=True)
    plt.plot([], [], color=adducts_color, linewidth=10,label="high-molecular-weight species", alpha=.3)
    plt.plot([], [], color=IgG_color, linewidth=10,label="native protein", alpha=.3)
    plt.xlim(ElutionRANGElow,ElutionRANGEhigh)
    plt.yticks([])
    # x-axis label
    # frequency label
    plt.ylabel('Absorbance')
    # plot title
    plt.title(seriesname)

#READING DATA IN
parser = argparse.ArgumentParser(description='File',prog='Analyze.py') #read the flags
parser.add_argument('-f', '--file', metavar='FILENAME.txt',dest='FILENAME',help='txt data file to be analyzed') # read flag for file into "FILENAME" variable
args = parser.parse_args()

if args.FILENAME==None:
    print "no file specified, analyzing files in ToProcess folder"
    os.listdir("./ToProcess")
    filelist = [f for f in os.listdir("./ToProcess") if os.path.isfile(os.path.join("./ToProcess", f))]
    filelist.sort(key=str.lower)
    plt.figure(figsize=(15,(3*len(filelist[:])))) # main plot
    for i in range(len(filelist[:])):
        file=filelist[i]
        subplot=(100*(len(filelist[:]))+10+i+1)
        #fig, ax = plt.subplots()
        ax=plt.subplot(subplot)
        (X,Y,name)=Test("./ToProcess/"+file)
        print "test done, now plotting"
        Plot(X,Y,name,colors[i])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.xaxis.set_ticks_position('bottom')
        if i==0:
            plt.legend(loc=2,fontsize='small')
        if i<len(filelist[:])-1: #not last plot
            ax.xaxis.set_major_formatter(plt.NullFormatter())
        else:
            plt.xlabel('Elution (mL)')
    
    plt.savefig("OUTPUT"+"_plot.png")
    plt.close()
    sys.exit()
elif not os.path.isfile(args.FILENAME):
    sys.exit(args.FILENAME+" does not exist. Please specify a DSC file, or leave blank and place file(s) to be analyzed in /ToProcess")
else:
    #Analyze(args.FILENAME)
    (X,Y,name)=Test(args.FILENAME)
    #plt.figure() # main plot
    fig, ax = plt.subplots()
    print "test done, now plotting"
    Plot(X,Y,name,colors[0])
    plt.savefig("EXAMPLE"+"_plot.png")
    plt.close()
    sys.exit()