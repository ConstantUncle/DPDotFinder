#import all necessary
import argparse
import subprocess
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import astropy.io.fits as fits
import scipy.interpolate
from astropy.table import Table, Column

#take arguments
parser = argparse.ArgumentParser()
parser.add_argument("filename", help="The event file or light curve to analyse")
parser.add_argument("period", help="The initial period to try analysing around", type=float)
parser.add_argument("pdiff", help="How much to increment the trial period by", type=float)
parser.add_argument("dpdot", help="The initial first derivative of the period to try analysing around", type=float)
parser.add_argument("dpdiff", help="How much to increment the trial first derivative by", type=float)
parser.add_argument("--trials", help="The number of iterations to try for each variable. Default is 7. Please note that this is iterated over twice, so the number of trials is the square of the number you enter.", default=7, type=int)
parser.add_argument("--cleanup", help="Whether DpdotFind should delete temporary files created during subroutines. Default is TRUE.", default=True, type=bool)
parser.add_argument("--outfile", help="Name of the output file to be created. Default is Results.", default="results")

args=parser.parse_args()

infile = args.filename
period = args.period
pdiff = args.pdiff
dpdot = args.dpdot
dpdiff = args.dpdiff
trials = args.trials
cleanup = args.cleanup
outfile = args.outfile

#Get info for output header
fil = fits.open(infile)
filehead = fil[1].header
TimeA = ("DATE-OBS", filehead["DATE-OBS"], filehead.comments["DATE-OBS"])
TimeB = ("TIME-OBS", filehead["TIME-OBS"], filehead.comments["TIME-OBS"])
TimeC = ("DATE-END", filehead["DATE-END"], filehead.comments["DATE-END"])
TimeD = ("TIME-END", filehead["TIME-END"], filehead.comments["TIME-END"])
RA = ("RA", filehead["RA_NOM"], filehead.comments["RA_NOM"])
DEC = ("DEC", filehead["DEC_NOM"], filehead.comments["DEC_NOM"])
fil.close()

#iterate efsearch
if trials%2==0:
    dpinit = dpdot - (trials/2 - 0.5)*dpdiff
    pinit = period - (trials/2 - 0.5)*pdiff
else:
    dpinit = dpdot - (trials//2)*dpdiff
    pinit = period - (trials//2)*pdiff

var = 0
err = 0
rows=[]
subprocess.call("mkdir DPDotTmp", shell=True)

for i in range(trials):
    dptry = dpinit + i*dpdiff
    for j in range(trials):
        periodtry = pinit + j*pdiff
        efs = "efsearch "+infile+" window='-' sepoch=INDEF dper="+str(periodtry)+" nphase=INDEF nbint=INDEF dres="+str(pdiff/128)+" nper=128 outfile=DPDotTmp/dpdot"+str(i)+str(j)+".fes plot=NO dpdot="+str(dptry)+" tchat=3"                                #command line string
        print("Iteration "+ str(i*trials + j+1)+" of "+str(trials**2))
        subprocess.call(efs, shell=True)                                              #running EFSearch
        t=Table.read("DPDotTmp/dpdot"+str(i)+str(j)+".fes", hdu=1, format="fits")     #Take table from file
        tfil = fits.open("DPDotTmp/dpdot"+str(i)+str(j)+".fes")                       #getting variance value for error
        var = tfil[1].header["VAROB_1"]
        err = np.sqrt(var)/128
        tfil.close()
        t.sort("CHISQRD1")                                                            #Sort by Chi2
        rows.append([t["PERIOD"][-1], dptry, t["CHISQRD1"][-1], err)                  #Append values to row

if cleanup==True:                                                                     #cleanup temp files
    subprocess.call("rm -r DPDotTmp", shell=True)

#Create results table and append rows
results = Table(rows=rows, names=["Period", "PDot", "Chi2", "Error"])

#output file
Hdu = fits.table_to_hdu(results)
head = Hdu.header
head.set("TUNIT1", "s", before="TTYPE2")
head.append(("INPUT", infile, "Input file path specified"), end=True)
head.append(("PIn", str(period), "Initial test period"), end=True)
head.append(("PDiff", str(pdiff), "Test Period step size"), end=True)
head.append(("DPDotIn", str(dpdot), "Initial test first derivative"), end=True)
head.append(("DPDiff", str(dpdiff), "Test first derivative step size"), end=True)
head.append(TimeA, end=True)
head.append(TimeB, end=True)
head.append(TimeC, end=True)
head.append(TimeD, end=True)
head.append(RA, end=True)
head.append(DEC, end=True)
Hdu.name = "Results"
Hdu.writeto(outfile+".dpd", overwrite=True)
#results.write("results.dpd", format="fits", overwrite=True)

#plot contours
results.sort("Chi2")
#plt.ion()

N=10000
x=(results["Period"])
y=(results["PDot"])
z=(results["Chi2"])

#tck = scipy.interpolate.bisplrep(x, y, z, s=0)

X = np.linspace(pinit, pinit+(trials-1)*pdiff, N)
Y = np.linspace(y.min(), y.max(), N)
#Z = scipy.interpolate.bisplev(X, Y, tck)
Z = scipy.interpolate.griddata((x, y), z, (X[None,:], Y[:,None]), method="cubic", rescale=True)

fig, ax = plt.subplots()
CS = ax.contourf (X, Y, Z, cmap="plasma")
#ax.clabel(CS)
ax.set_xlabel("Period", fontsize=8)
ax.set_ylabel("PDot", fontsize=8)
ax.set_title("Best Period="+str(results["Period"][-1])+", Best PDot="+str(results["PDot"][-1]), fontsize=10, pad=12)
fig.colorbar(CS)
plt.show()

#print best values
print("Best Period="+str(results["Period"][-1])+", Best PDot="+str(results["PDot"][-1]))
print("Error="+str(results["Error"][-1]))

