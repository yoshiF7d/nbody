import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
import struct
import subprocess
import matplotlib.style as mplstyle

mplstyle.use('fast')
EPSILON = 0.1

def parselim(limstr):
	l0,l1 = limstr.split(",")
	l0 = l0.lstrip('([').lstrip()
	l1 = l1.rstrip(')]').rstrip()
	return [float(l0),float(l1)]
	
def kin(m,px,py):
	return np.sum(0.5*(np.square(px) + np.square(py)) / m)
	
def pot(c,x,y):
	energy=0
	for i in range(c.size):
		r = np.delete(np.sqrt(np.square(x-x[i]) + np.square(y-y[i]) + EPSILON),i)
		c2 = np.delete(c[i]*c,i)
		energy += np.sum(c2/r)
			
	return energy

parser = argparse.ArgumentParser(description='plot nbody simulation result')
parser.add_argument('indir')
parser.add_argument('outdir',default='.',nargs='?')
parser.add_argument('mp4file',nargs='?')
parser.add_argument('enfile',default='energy.txt',nargs='?')
parser.add_argument('-f','--flim',type=str,help='process from fmin th file to fmax th file (index starts from 0). usage : -f [fmin,fmax]',metavar='[fmin,fmax]')
args=parser.parse_args()

if args.flim is not None:
	args.flim = list(map(int,parselim(args.flim)))
	files = files[args.flim[0]:args.flim[1]]

files=[file for file in os.listdir(args.indir) if file.endswith('.data')]
files.sort()

data=[]

for file in files:
	dtype=[('time','f8'),('id','i4'),('mass','f8'),('charge','f8'),('x','f8'),('y','f8'),('px','f8'),('py','f8'),('fx','f8'),('fy','f8')]
	data.append(np.fromfile(os.path.join(args.indir,file), dtype=dtype))

t = np.empty(len(data))
ken = np.empty(len(data))
pen = np.empty(len(data))

for i,d in enumerate(data):
	t[i] = d['time'][0]
	ken[i] = kin(d['mass'],d['px'],d['py'])
	pen[i] = pot(d['charge'],d['mass'],d['x'])

#s = struct.Struct('=di8d')
#with open(os.path.join(args.indir,files[0]),'rb') as f:
#	data = f.read()

if not os.path.exists(args.outdir):
	os.makedirs(args.outdir)

fig,axs=plt.subplots(1,2,figsize=(12, 6))
#particles plot
axs[0].set_xlim(-10,10)
axs[0].set_ylim(-10,10)
axs[0].set_aspect('equal')
#energy plot
ymax=np.max(ken+pen)
xmax=t[-1]
axs[1].set_xlim(0,xmax)
axs[1].set_ylim(0,ymax)
axs[1].set_aspect(xmax/ymax)
lines=[None]*3

for i,file in enumerate(files):
	infile = os.path.join(args.indir,file)
	outfile = os.path.join(args.outdir,os.path.splitext(os.path.basename(file))[0]+'.png')
	
	print(infile + ' -> ' + outfile)
	
	if i==0:
		scat=axs[0].scatter(data[i]['x'],data[i]['y'],c=data[i]['charge'],cmap='bwr')
		lines[0],=axs[1].plot(t[:i],ken[:i],label='kinetic')
		lines[1],=axs[1].plot(t[:i],pen[:i],label='potential')
		lines[2],=axs[1].plot(t[:i],ken[:i]+pen[:i],label='total')
		axs[1].legend()
	else:
		scat.set_offsets(list(zip(data[i]['x'],data[i]['y'])))
		lines[0].set_data(t[:i],ken[:i])
		lines[1].set_data(t[:i],pen[:i])
		lines[2].set_data(t[:i],ken[:i]+pen[:i])
	
	fig.suptitle('time : {:>.2e}'.format(t[i]))
	fig.savefig(outfile,format='png',dpi=80)

with open(args.enfile,'w') as f:
	f.write('time\tken\tpen\tten\n')
	for i in range(len(data)):
		f.write(f'{t[i]}\t{ken[i]}\t{pen[i]}\t{ken[i]+pen[i]}\n')

if args.mp4file:
	args = ['ffmpeg','-r','30','-i',os.path.join(args.outdir,'%04d.png'),'-vcodec','libx264','-pix_fmt','yuv420p','-r','30',args.mp4file]
	subprocess.run(args)

