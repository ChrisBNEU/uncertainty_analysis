from mpi4py import MPI
import time
import sys

def test_mpi4py():
    start = time.asctime()
    size = MPI.COMM_WORLD.Get_size()
    rank = MPI.COMM_WORLD.Get_rank()
    name = MPI.Get_processor_name()
    
    time.sleep(5)
    stop = time.asctime()
    
    sys.stdout.write(
        "Hello, World! I am process %d of %d on %s \n started at %s\n, finished at %s.\n"
        % (rank, size, name, start, stop))
    
if (__name__ == "__main__"):
    test_mpi4py()
    