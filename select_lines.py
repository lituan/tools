import os
import sys

def read_id(id_f):
    with open(id_f) as o_f:
        lines = o_f.readlines()
        lines = [line.rstrip('\r\n') for line in lines]
        lines = [line for line in lines if line]

        return lines

def read_check(check_f,ids):
    with open(check_f) as o_f:
        lines = o_f.readlines()
        lines = [line.rstrip('\r\n') for line in lines]
        lines = [line for line in lines if line]
        check_lines = []
        for line in lines:
            for i in ids:
                if i in line:
                    check_lines.append(line)
    return check_lines

def main():
    ids = read_id(sys.argv[-2])
    check_lines = read_check(sys.argv[-1],ids)
    with open('check.txt','w') as w_f:
        for line in check_lines:
            print >> w_f,line

if __name__ == "__main__":
    main()


