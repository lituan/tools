#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
input command and directory containing data files, these data files are read at the beginning
and not read a second time, then this script will find available queue and submit jobs

default walltime is 240 hours

example
python auot_qsub.py 'wc' temp
apply command 'wc' to all files contained in temp, and qsub every single jobs
"""
import sys
import os
import subprocess

def queue_state():
    lines = subprocess.check_output(['statj'], shell=True)
    lines = lines.splitlines()

    # check available state
    good_queue = []
    for line in lines:
        if 'available cores' in line:
            words = line.split()
            available = float(words[-1])
            good_queue.append((words[0],available))
    # check waiting state
    waiting = []
    for line in lines:
        if 'total' in line:
            words = line.split()
            if len(words) == 3:
                waiting.append(float(words[-1]))

    #find good queue, no waiting and more available
    good_queue = [(q[1]-waiting[i],q[0]) for i,q in enumerate(good_queue)]
    good_queue = sorted(good_queue,reverse=True)
    good_queue = [q[-1] for q in good_queue]

    # ignore q5, specially for GPU
    if 'q5' in good_queue:
        good_queue.pop(good_queue.index('q5'))

    return good_queue


def files_in_dir(directory):
    final_files = []
    for root, dirs, files in os.walk(directory):
        for f in files:
            final_files.append(os.path.join(root, f))
    return final_files


def create_temp_pbs(command, queue,temp_name):
    with open(temp_name, 'w') as w_f:
        print >> w_f, '#! /bin/bash'
        print >> w_f, '#PBS -q' + ' ' + queue
        print >> w_f, '#PBS -l nodes=1:ppn=1,walltime=240:00:00' # hour:minute:seconds
        print >> w_f, '#PBS -j oe'
        print >> w_f, 'cd $PBS_O_WORKDIR'
        print >> w_f, ' '
        print >> w_f, '#RUN'
        print >> w_f, command


def delete_temp_pbs(temp_name):
    subprocess.call('rm '+temp_name, shell=True)


def main():
    command = sys.argv[1]
    total_command = ' '.join(sys.argv[1:])
    f = sys.argv[-1]
    f_name = os.path.splitext(os.path.split(f)[-1])[0]
    log_f = open(f_name + "_" + 'auto_qsub_log.txt', 'w')
    temp_name = command+'_'+f_name
    # try to submmit the job to the available queue by trying each queue by
    # the order
    while len(good_queue) < 1:
        good_queue = queue_state()
    while len(good_queue) > 0:
        queue = good_queue[0]
        create_temp_pbs(total_command, queue,temp_name)
        try:
            output = subprocess.check_output(['qsub', temp_name])
            lines = output.splitlines()
            output = '\n'.join(lines)
            print 'qsub success',queue
            print total_command
            print output
            print >> log_f,queue
            print >> log_f, output
            print >> log_f, total_command
            print >> log_f, '*' * 80
            delete_temp_pbs(temp_name)
            break # if success, stop trying rest queues and try next data_file
        except Exception,e:
            print e
            continue
        finally:
            good_queue = queue_state()
            while len(good_queue) < 1:
                good_queue = queue_state()

    log_f.close()

if __name__ == "__main__":
    main()


