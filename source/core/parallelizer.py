import multiprocessing as mp
from progressbar import ProgressBar
from collections import OrderedDict

class Parallelizer(object):

    def __init__(self, workers=-1):
        self.workers = mp.cpu_count() if workers == -1 else workers
        self.queue = mp.Queue()
        self.results = None

    def prun(self, function, data):
        self.N = len(data)
        self.tasks = []
        chunksize = 1.*self.N / self.workers
        try:
            for i in xrange(self.workers):
                start = int(round(i*chunksize))
                end = int(round((i+1)*chunksize))
                datachunk = data[start:end]
                p = mp.Process(target = function,
                                 args = (datachunk, start, self.queue, ))
                self.tasks.append(p)
                p.start()
            accum = {}
            with ProgressBar(self.N) as pb:
                for j in xrange(self.N):
                    result = self.queue.get(block=True, timeout=None)
                    accum.update(result)
                    pb.update(j)
            ordered_results = OrderedDict(sorted(accum.items()))
            self.results = ordered_results.values()
        except (KeyboardInterrupt, SystemExit):
            for task in self.tasks:
                task.terminate()
        except Exception, e:
            print "**Caught exception: ", e
            print "\nTerminating workers...\n"
            for task in self.tasks:
                task.terminate()
        finally:
            for task in self.tasks:
                task.join()
            self.accumulator.close()