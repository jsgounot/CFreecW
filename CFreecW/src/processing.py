import os, tempfile, time, sys
from multiprocessing import Pool, cpu_count
from subprocess import check_output

from Bio import SeqIO

class TMPFname() :

    def __init__(self, delete=True, ext="", ** kwargs) :
        suffix = self.format_ext(ext)
        tmp = tempfile.NamedTemporaryFile(delete=False, suffix=suffix, ** kwargs)
        self.fname = tmp.name 
        self.delete = delete

    def __str__(self) :
        return self.fname

    def format_ext(self, ext) :
        ext = str(ext)
        if not ext : return ext
        if not ext.startswith(".") : ext = "." + ext.strip()
        return ext

    def exist(self) :
        return os.path.isfile(self.fname)

    def remove(self) :
        if self.exist() and self.delete :
            print ("Delete temporary file at : %s" %(self.fname))
            os.remove(self.fname)

    def __del__(self) :
        self.remove()

class TMPFasta(TMPFname) :

    def __init__(self, data, delete=True, wodesc=False, ** kwargs) :
        super(TMPFasta, self).__init__(delete=delete, ext="fa", ** kwargs)
        if data : self.write(data, wodesc=wodesc)

    @staticmethod
    def remove_desc(feature) :
        feature.description = ""
        return feature

    def write(self, data, mode="w", wodesc=False) :
        if wodesc : data = (TMPFasta.remove_desc(feat) for feat in data)
        with open(self.fname, mode) as handle :
            SeqIO.write(data, handle, "fasta")


"""
PROCESSING
"""

# decorator to time functions
def ftiming(f):

    def wrap(* args, ** kwargs):      
        
        funcname = kwargs.pop("funcname", f.__name__)

        time1 = time.time()
        ret = f(* args, ** kwargs)
        time2 = time.time()
        
        print('{:s} function took {:.3f} ms'.format(funcname, (time2-time1) * 1000.0))
        
        return ret

    return wrap

@ftiming
def run_cmdline(cmdline, funcname=None) :
    return check_output(cmdline, universal_newlines=True, shell=True)

def readfile(fname) :
    if not os.path.isfile(fname) :
        raise OSError("File not found : %s" %(fname))
    
    with open(fname) as f :
        for line in f :
            line = line.strip()
            yield line

"""
MULTIPROCESSING
"""

class FunArgs() :
 
    def __init__(self, fun, * args, ** kwargs) :
        self.fun = fun
        self.args = args
        self.kwargs = kwargs
 
    def launch_fun(self) :
        return self.fun(* self.args, ** self.kwargs)
 
class Multiprocess() :
 
    @staticmethod
    def lambda_fun(farg) :
        return farg.launch_fun()
 
    def run(self, fargs, ncore=1) :
        print ("ncore : %i - available : %i" %(ncore, cpu_count()))
        if ncore == 1 : return [farg.launch_fun() for farg in fargs]
 
        pool = Pool(ncore)
        func = Multiprocess.lambda_fun
        fargs = ((farg, ) for farg in fargs)
 
        try :
            data = pool.starmap(func, fargs)
            pool.close()
            pool.join()
            return data
 
        except KeyboardInterrupt :
            print("Interrupt childs process")
            pool.terminate()
            sys.exit()
 
def mproc(fargs, ncore=1) :
    mp = Multiprocess()
    return mp.run(fargs, ncore=ncore)