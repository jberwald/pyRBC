import multiprocessing as mp
import time
import numpy as np
import itertools

def foo_pool(x):
    time.sleep(2)
    print 'called'
    return x*x

result_list = []
def log_result(result):
    # This is called whenever foo_pool(i) returns a result.
    # result_list is modified only by the main process, not the pool workers.
    result_list.append(result)

def apply_async_with_callback(f):
    #input a function (ex foo_pool)
    pool = mp.Pool()
    for i in range(10):
        pool.apply_async(f, args = (i, ), callback = log_result)
    pool.close()
    pool.join()
    print(result_list)

def test_apply_async(f,arrs):
    result_list = []
    pool = mp.Pool()
    list = arange(0,10)
    pool.map (f, list)
    print(result_list)