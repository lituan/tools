import sys
import os
import lt

def write_pbs():
    a = '''#! /bin/bash
#PBS -q q1
#PBS -l nodes=1:ppn=1,walltime=24:00:00
#PBS -j oe
cd $PBS_O_WORKDIR

#RUN
'''
    for i,f in enumerate(lt.files_in_dir(sys.argv[-1])):
        ofile = lt.openfile('pbs_'+str(i),fileextention='.pbs')
        print >> ofile,a
        print >> ofile,'{0}\t~/{1}\t~/{2}\t~/{3}'.format('python',sys.argv[-2],f,'wd648.wdsp')

write_pbs()


