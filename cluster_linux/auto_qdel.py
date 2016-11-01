#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
auto delete
"""

import os
import sys
import subprocess

queues = subprocess.check_output('qstat -u lt')
queues = queues.splitlines()
queues = [q.split()[0].split('.')[0] for q in queues if q.split()[2] == 'lt']

for q in queues:
    command = 'qdel' + ' ' + q
    subprocess.check_output(command,shell=True)
