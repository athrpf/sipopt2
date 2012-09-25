import subprocess
import sys
import random

class Interval:    # 1 intervall: upper and lower bound
    def __init__(self, pL, pU): # standard constructor
        assert(len(pL)==len(pU))
        self.pL = pL
        self.pU = pU

    def split(self, idx): # split an interval -- in: 1 interval, out: 2 intervals resulting from split
        p1L = []
        p1U = []
        p2L = []
        p2U = []
        for k, (pLv, pUv) in enumerate(zip(self.pL, self.pU)):
            p1L.append(pLv)
            p2U.append(pUv)
            if k==idx:
                v = 0.5*(pLv+pUv)
                p1U.append(v)
                p2L.append(v)
            else:
                p1U.append(pUv)
                p2L.append(pLv)
        return Interval(p1L, p1U), Interval(p2L, p2U)

class AmplSet: # an ampl setup: ampl srcipt and feasible set of intervals
    def __init__(self, ampl_script, pLnames, pUnames, pLvalues, pUvalues):
        self.ampl_script = ampl_script
        self.pLnames = pLnames
        self.pUnames = pUnames
        self.pLvalues = pLvalues
        self.pUvalues = pUvalues
        self.intervals = [Interval(pLvalues, pUvalues)]

    def write_include_file(self, filename='intervalls.inc'): # write current intervallization data into ampl include file
        f = open(filename, 'w')
        nint = len(self.intervals)
        f.write('let nint := {};\n'.format(nint))
        for name in self.pLnames:
            f.write('let {{q in 1..{} }} {}[q].intervalID := q;\n'.format(nint, name))
        for k, iv in enumerate(self.intervals):
            for nameL, pL, nameU, pU in zip(self.pLnames, iv.pL, self.pUnames, iv.pU):
                f.write('let {}[{}] := {};\n'.format(nameL, k+1, pL))
                f.write('let {}[{}] := {};\n'.format(nameU, k+1, pU))
        f.close()

    def call_ampl(self): # one run of ampl: write inc file and run ampl with the specified setup
        self.write_include_file()
        try:
            subprocess.check_call(['ampl', self.ampl_script])
        except subprocess.CalledProcessError as e:
            print e
            sys.exit(1)

    def split(self, interval_idx, para_idx): # split interval at interval_idx with respect to parameter para_idx
        split_interval = self.intervals.pop(interval_idx)
        interval1, interval2 = split_interval.split(para_idx)
        self.intervals.extend([interval1, interval2])

    def randomize(self,nruns):
        for ir in range(nruns):
            self.call_ampl()
            nint = random.randint(0,(len(self.intervals)-1))
            npar = random.randint(0,(len(self.pLnames)-1))
            self.split(nint,npar)


def run():
    ampl_script = 'autorandomintervals.run'
    pLnames = ['p1L', 'p2L']
    pUnames = ['p1U', 'p2U']
    pLvalues = [0.9, 0.8]
    pUvalues = [1.1, 1.2]
    info = AmplSet(ampl_script, pLnames, pUnames, pLvalues, pUvalues)
    info.randomize(3)

if __name__=='__main__':
    run()
