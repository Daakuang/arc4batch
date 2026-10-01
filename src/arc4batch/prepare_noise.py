"""Create or verify the exact paired noise streams used in the study."""
import argparse
import numpy as np
from scipy.io import savemat,loadmat
from adaptive_model import OUT

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--seeds',type=int,nargs='+',default=[0,50,90,94,191,200,201,202,203,204]);a=p.parse_args()
    OUT.mkdir(parents=True,exist_ok=True)
    for seed in a.seeds:
        expected=np.random.Generator(np.random.PCG64(seed)).standard_normal((5,30001))
        path=OUT/f'noise_{seed:03d}.mat'
        if path.exists():assert np.array_equal(loadmat(path)['paired_noise'],expected),str(path)
        else:savemat(path,dict(paired_noise=expected))
    print('Noise streams created or verified:',a.seeds)
