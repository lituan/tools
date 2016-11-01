import datetime
import itertools
import operator
import os
import random
import re
import sys
import time
from functools import wraps
from collections import OrderedDict
import cPickle as pickle


def zen_of_python():

    ZEN = ['Beautiful is better than ugly',
           'Explicit is better than implicit',
           'Simple is better than complex',
           'Complex is better than complicated',
           'Flat is better than nested',
           'Sparse is better than dense',
           'Readability counts',
           "Special cases aren't special enough to break the rules",
           'Although practicality beats purity',
           'Errors should never pass silently',
           'Unless explicitly silenced',
           'In the face of ambiguity, refuse the temptation to guess',
           'There should be one-- and preferably only one -- obvious way to do it',
           "Although that way may not be obvious at first unless you're Dutch",
           'Now is better than never',
           'Although never is often better than *right* now',
           "If the implementation is hard to explain, it's a bad idea",
           'If the implementation is easy to explain, it may be a good idea',
           "Namespaces are one honking great idea -- let's do more of those!"]

    print '-' * 80
    print datetime.datetime.now()
    print ZEN[random.randint(0, 18)]
    print '-' * 80


def log(fn):
    @wraps(fn)
    def wrap(*args, **kwds):
        with cd.open(fn.__name__ + '.log', 'w', encoding='utf_8') as log_f:
            print >> log_f, '# -*- coding:utf-8 -*-'
            back_stderr, back_stdout = sys.stderr, sys.stdout
            sys.stderr, sys.stdout = log_f, log_f
            result = fn(*args, **kwds)
            sys.stderr, sys.stdout = back_stderr, back_stdout
            return result
    return wrap


def run_time(fn):
    @wraps(fn)
    def wrapper(*arg, **kwds):
        start = time.time()
        try:
            return fn(*arg, **kwds)
        finally:
            zen_of_python()
            end = time.time()
            process = (end - start)
            if process > 100:
                process = process / 60.0
                print '{0:s} run {1:.2f}minutes'.format(fn.__name__, process)
            elif process < 1:
                process = process * 1000
                print '{0:s} run {1:.2f}us'.format(fn.__name__, process)
            else:
                print '{0:s} run {1:.2f}s'.format(fn.__name__, process)
            print '-' * 80
    return wrapper


def print_list(lis, n=1, output=sys.stdout):
    print '*' * 80
    for l in [lis[i:i + n] for i in range(0, len(lis), n)]:
        print ' '.join(l)


def print_dic(dic, nest=0, output=sys.stdout):
    spacing = '    '*2
    for k, v in dic.iteritems():
        if type(v) == dict or type(v) == OrderedDict:
            print >> output, '{0}{1}'.format((nest) * spacing, k)
            print_dic(v, nest + 1, output)
        elif type(v) == list:
            print >> output, '{0}{1}'.format((nest) * spacing, k)
            for i in v:
                print >> output, '{0}{1}'.format((nest + 1) * spacing, i)
        else:
            print >> output, '{0}{1}'.format((nest) * spacing, k)
            print >> output, '{0}{1}'.format((nest + 1) * spacing, v)


def open_file(file_name='', file_suffix='', file_extension='.txt', dir_name='', dir_suffix='', inner_dir='', result_path='', mod='w', argv=sys.argv):

    script_path, script_name = os.path.split(argv[0])
    script_name, script_extension = os.path.splitext(script_name)
    if len(argv) == 1:
        f_path, f_name = '', ''
    if len(argv) == 2:
        f_path, f_name = os.path.split(argv[-1])
        f_name, _ = os.path.splitext(f_name)
    elif len(argv) >= 3:
        if ( os.path.isfile(argv[-1]) or os.path.isdir(argv[-1]) ) and ( os.path.isfile(argv[-1]) or os.path.isdir(argv[-2]) ):
            f_path1, f_name1 = os.path.split(argv[-2])
            f_name1, _ = os.path.splitext(f_name1)
            f_path2, f_name2 = os.path.split(argv[-1])
            f_name2, _ = os.path.splitext(f_name2)
            f_path, f_name = f_path1, f_name1 + '_' + f_name2
        else:
            f_path, f_name = os.path.split(argv[-1])
            f_name, _ = os.path.splitext(f_name)

    if file_name == '':
        if file_suffix == '':
            if f_name == '':
                result_name = script_name
            else:
                result_name = f_name + '_' + script_name
        else:
            if f_name == '':
                result_name = script_name + '_' + file_suffix
            else:
                result_name = f_name + '_' + script_name + '_' + file_suffix
    else:
        if file_suffix == '':
            result_name = file_name
        else:
            result_name = file_name + '_' + file_suffix

    if dir_name == '':
        if dir_suffix == '':
            if f_name == '':
                result_dir = script_name
            else:
                result_dir = f_name + '_' + script_name
        else:
            if f_name == '':
                result_dir = script_name + '_' + dir_suffix
            else:
                result_dir = f_name + '_' + script_name + '_' + dir_suffix
    else:
        if dir_suffix == '':
            result_dir = dir_name
        else:
            result_dir = dir_name + '_' + dir_suffix

    if result_path == '':
        if f_path == '':
            result_path = script_path
        else:
            result_path = f_path
    else:
        result_path = reslut_path

    result_path = os.path.join(result_path, result_dir)

    if inner_dir != '':
        result_path = os.path.join(result_path, inner_dir)

    result_path = re.sub(u'[*?"<>|]', '-', result_path)
    if not os.path.exists(result_path):
        os.makedirs(result_path)


    result_file = os.path.join(result_path, result_name + file_extension)
    result_file = re.sub(u'[*?"<>|]', '-', result_file)
    if not os.path.exists(result_file):
        file_handle = open(result_file, mod)
        return file_handle
    else:
        # replace existed file
        file_handle = open(result_file, mod)
        return file_handle
        for i in range(1000000):
            result_file = result_name + '_v' + str(i + 1)
            result_file = os.path.join(
                result_path, result_file + file_extension)
            result_file = re.sub(u'[:*?"<>|]', '-', result_file)
            if not os.path.exists(result_file):
                file_handle = open(result_file, mod)
                return file_handle
            else:
                continue


def files_in_dir(directory):
    for root, dirs, files in os.walk(directory):
        for f in files:
            yield os.path.join(root, f)


def align_lis_lis(lis_lis):
    """align and trans nested list to print a table"""
    lis_lis = [[str(l) for l in lis]
               for lis in lis_lis]  # trans every element to str
    #make all inner lists of the same length
    inner_lis_max_len = max(len(lis) for lis in lis_lis)
    lis_lis = [lis + (inner_lis_max_len - len(lis)) * [''] for lis in lis_lis]
    #trans list, so that the elements of the same column are in one list
    lis_lis = [[lis[i] for lis in lis_lis] for i in range(inner_lis_max_len)]
    #make element in the same list have the same length
    aligned = []
    for lis in lis_lis:
        width = max([len(l) for l in lis])
        lis = [l + (width - len(l)) * ' ' for l in lis]
        aligned.append(lis)
    #trans list_list to the original list_list
    inner_lis_max_len = max(len(lis) for lis in lis_lis)
    lis_lis = [[lis[i] for lis in aligned] for i in range(inner_lis_max_len)]
    return lis_lis

# a = [['dd','dddd'],['ddddd','dd']]
# print align_lis_lis(a)
# a = [['dd  ','dddd'],['ddddd','dd  ']]


def trans_dic_lis(dic_lis):
    lis_lis = [[k] + v for k, v in dic_lis.iteritems()]
    return trans_lis_lis(lis_lis)


def lis_sta(lis):
    key = set(lis)
    sta = [(i, lis.count(i)) for i in key]
    sta = sorted(sta, reverse=True)
    return sta


def read(f):
    with open(f) as o_f:
        lines = o_f.readlines()
        lines = [line for line in lines if len(line.split()) > 0]
        lines = [line.strip('\n') for line in lines]
        return lines

def pickle_dump(obj,f_name):
    pickle.dump(obj,open_file(file_name=f_name,file_extension='.pickle'))

def pickle_load(obj_f):
    return pickle.load(open(obj_f,'r'))

