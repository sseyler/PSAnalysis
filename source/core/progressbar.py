'''
Created on Nov 20, 2012

@author: Sean Lee Seyler
'''
import time
import sys
import numpy as np

class ProgressBar(object):
    def __init__(self, N, width=50, update_interval=1, unit='secs', verbose=True):
        # initialize global features
        self.N = N
        self.status = ""
        # initialize bar features
        self.bar_width = width
        # initialize eta features
        self.t_out_interval = update_interval
        self.etatxt = "ETA: ..."
        # initialize wall time features
        self.unit = unit
        self.verbose = verbose

    def update(self, progress):
        if isinstance(progress, int) or isinstance(progress, float):
            fraction = float(progress)/self.N
        else:
            fraction = 0
            self.status = "Error: progress variable must be float or int\r\n"
        if fraction < 0:
            fraction = 0
            self.status = "Halt...\r\n"
        if fraction >= 1:
            fraction = 1
            self.status = "Done...\r\n"
        bar = self._update_text(fraction, time.time())
        sys.stdout.write(bar)
        sys.stdout.flush()

    # Find the best estimator for remaining time. Simplest solution is to use the weighted
    # average of:
    #   1) the average overall completion rate
    #   2) the average completion rate over some moving interval (window), e.g., the
    #       the completion rate over the last x seconds.
    # Progress variable rescaled automatically w/ self.N (length of system to be monitored)
    def _update_text(self, fraction, t):
        if t >= self.t_out:
            try:
                frac_left = 1 - fraction
                t_elap = t - self.start
                avg_rate = t_elap/fraction
                t_left = np.rint(frac_left*avg_rate)
                mins = t_left/60
                secs = t_left%60
                self.etatxt = "ETA: %im%is" % (mins, secs)
            except ZeroDivisionError:
                self.etatxt = "ETA: ..."
            finally:
                self.t_out += self.t_out_interval
                self.frac_intvrl_1 = fraction
        block = np.rint(self.bar_width*fraction)
        return "\r[%s] %4.1f%% | %s %s" % ( "#"*block + "-"*(self.bar_width-block), fraction*100, self.etatxt, self.status)     

    def __enter__(self):
        # setup toolbar
        sys.stdout.write("\nProgress:\n")
        sys.stdout.write("[%s]" % ("-" * self.bar_width))
        sys.stdout.flush()
        sys.stdout.write("\b" * (self.bar_width+1)) # return to start of line, after '['
        sys.stdout.flush()

        # start timing
        self.start = time.time()
        self.t_out = self.start + self.t_out_interval
        # interval 1 is beginning of update interval, interval 2 is 3x longer than interval 1 or
        # 5 seconds, whichever is longers
        # self.t_intrvl = self.start
        self.frac_intvrl = 0
        return self

    def __exit__(self, *args):
        self.end = time.time()

        sys.stdout.write("\n") # Add newline after progress bar
        sys.stdout.flush() # Reset line

        self.secs = self.end - self.start
        if self.verbose:
            if self.unit == 'msecs':
                self.msecs = self.secs * 1000  # millisecs
                print '  --->  Elapsed time: %.0f ms\n' % self.msecs
            elif self.unit == 'secs':
                print '  --->  Elapsed time: %.2f s\n' % self.secs