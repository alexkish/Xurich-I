from ROOT import *
from math import sqrt, pi, fabs
from array import array
from operator import add
import sys

#import psyco
#psyco.full()

#import pdb
#pdb.set_trace()

#set the detector dimensions
sensitive_volume_radius = 100.0 # mm
sensitive_volume_length = 100.0 # mm
Zoffset = 50.0

#DATA PATH
#schroedinger
#root_files_directory = '/home/physik/alexkish/data/arx/'
#xecluster
#root_files_directory='/archive/mc/xenon100/alexkish/CalibrationSources/'
#root_files_directory='/archive/mc/xenon100/alexkish/bg/t1/LXe/' 
#idark05x
root_files_directory='/Users/alexkish/data/xuerich/mc/compton/'

root_file = sys.argv[1]

filename = root_files_directory + root_file
file = TFile(filename, 'read')
tree = file.Get('t1')

nb_entries = tree.GetEntries()
#nb_entries = 10000

newfilename = root_files_directory
newfile = TFile(root_files_directory + root_file[:-5] + '-t2_0.1mm.root', 'recreate')
newtree = TTree('t2', '')

#che eto takoe??
max_size = 100

eventid	= array('i', [0])
etot	= array('f', [0.])
ns		= array('i', [0])
ed 		= array('f', max_size*[0.])
xp 		= array('f', max_size*[0.])
yp 		= array('f', max_size*[0.])
zp 		= array('f', max_size*[0.])
etotv 	= array('f', [0.])
nsv 	= array('i', [0])
edv 	= array('f', max_size*[0.])
xpv 	= array('f', max_size*[0.])
ypv 	= array('f', max_size*[0.])
zpv 	= array('f', max_size*[0.])
xp_pri 	= array('f', [0.])
yp_pri 	= array('f', [0.])
zp_pri 	= array('f', [0.])

newtree.Branch('eventid', 	eventid,	'eventid/I')
newtree.Branch('etot', 		etot, 		'etot/F')
newtree.Branch('ns', 		ns, 		'ns/I')
newtree.Branch('ed', 		ed, 		'ed[ns]/F')
newtree.Branch('xp', 		xp, 		'xp[ns]/F')
newtree.Branch('yp', 		yp, 		'yp[ns]/F')
newtree.Branch('zp', 		zp, 		'zp[ns]/F')
newtree.Branch('xp_pri', 	xp_pri, 	'xp_pri/F')
newtree.Branch('yp_pri', 	yp_pri, 	'yp_pri/F')
newtree.Branch('zp_pri', 	zp_pri, 	'zp_pri/F')

for entry in range(nb_entries):
#for entry in range(1000):
	tree.GetEntry(entry)

	if entry % 100 == 0:
	#if entry % 10 == 0:
		print '(%d/%d)' % (entry, nb_entries)
		#print '(%d/%d)' % (entry, 1000)
	
	tree_xp, tree_yp, tree_zp = list(tree.xp), list(tree.yp), list(tree.zp)
	tree_ed = list(tree.ed)
	tree_type = list(tree.type)

	# build a list of energy depositions [[e0, [x0, y0, z0]], ...] for all steps, sorted in z
	all_steps = [[tree_ed[step], [tree_xp[step], tree_yp[step], tree_zp[step]], tree_type[step]] for step in range(tree.nsteps)]
	all_steps.sort(lambda p1, p2: int(p1[1][2]-p2[1][2]))

	# select only steps in the sensitive volume and with energy deposited
	filter_er = lambda p: p[2] == 'gamma' or p[2] == 'e-' or p[2] == 'e+'
    #filter_er = lambda p: p[2] == 'e+'
	filter_ed = lambda p: p[0] > 0.
	filter_zp = lambda p: p[1][2] < Zoffset and p[1][2] > -sensitive_volume_length+Zoffset
	filter_rp = lambda p: sqrt(p[1][0]**2 + p[1][1]**2) < sensitive_volume_radius
	er_steps = filter(filter_er, all_steps)
	nz_steps = filter(filter_ed, er_steps)
	sv_steps = filter(filter_zp, nz_steps)
	sv_steps = filter(filter_rp, sv_steps)


	# compute the number of scatters in the sensitive volume
	if len(sv_steps) > 0:
		# start with 1 scatter
		nb_s = 1

		s_etot = sv_steps[0][0]
			
		s_ed = sv_steps[0][0]
		s_xp = sv_steps[0][1][0]*sv_steps[0][0]
		s_yp = sv_steps[0][1][1]*sv_steps[0][0]
		s_zp = sv_steps[0][1][2]*sv_steps[0][0]

		# for all steps starting with the second one
		for step in range(1, len(sv_steps)):
			step_ed = sv_steps[step][0]
			step_xp, step_yp, step_zp = sv_steps[step][1][0], sv_steps[step][1][1], sv_steps[step][1][2]

			# if the step if separated enough
			#if fabs(step_zp - sv_steps[step-1][1][2]) > 3.0:
			if fabs(sqrt((step_xp - sv_steps[step-1][1][0])**2 + (step_yp - sv_steps[step-1][1][1])**2 + (step_zp - sv_steps[step-1][1][2])**2)) > 0.1:
				# then this is a new scatter

				# save the energy deposited and the weighted position of the previous scatter
				ed[nb_s-1] = s_ed
				xp[nb_s-1] = s_xp/s_ed
				yp[nb_s-1] = s_yp/s_ed
				zp[nb_s-1] = s_zp/s_ed

				s_etot += step_ed

				# start a new scatter
				nb_s += 1
				s_ed = step_ed
				s_xp, s_yp, s_zp = step_xp*step_ed, step_yp*step_ed, step_zp*step_ed

			else:
				# this is not a new scatter, keep adding in the current
				s_etot += step_ed
				s_ed += step_ed
				s_xp += step_xp*step_ed
				s_yp += step_yp*step_ed
				s_zp += step_zp*step_ed

		# save the last scatter (or the first if we have only one!)
		ed[nb_s-1] = s_ed
		xp[nb_s-1] = s_xp/s_ed
		yp[nb_s-1] = s_yp/s_ed
		zp[nb_s-1] = s_zp/s_ed
	else:
		nb_s = 0


	if nb_s>0:
		# save event information
		eventid[0] = tree.eventid

		etot[0] = s_etot
		ns[0] = nb_s

		xp_pri[0] = tree.xp_pri
		yp_pri[0] = tree.yp_pri
		zp_pri[0] = tree.zp_pri

		# fill the tree
		newtree.Fill()

	else:
		eventid[0] = tree.eventid
		etot[0] = 0
		ns[0] = 0

		newtree.Fill()

file.Close()

newfile.Write()
newfile.Close()

