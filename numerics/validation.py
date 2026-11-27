from pathlib import Path
import csv
import numpy as np
import verify as b
from study import fit

def exact(v,w,margin):
    D,loss,L=b.D,b.LOSS,b.L
    x=np.array([L,2*L])[:,None]
    s=np.array([0.,3*L])[None,:]
    k=np.sqrt((v/(2*D))**2+(loss+1j*w)/D)
    extent=2*margin+3*L
    return (np.exp(v*(x-s)/(2*D))*np.sinh(k*(np.minimum(x,s)+margin))
            *np.sinh(k*(3*L+margin-np.maximum(x,s)))/(D*k*np.sinh(k*extent)))

def main():
    rows=[]
    for h in [.1,.05,.025]:
        errors=[]
        for v in [.1,-.1,0.]:
            for w in b.W:
                target=exact(v,w,20.)
                errors.extend((abs(b.forward(w,v,20.,h)/target-1)).ravel())
        rows.append((h,max(errors)))
    with Path(__file__).with_name('grid.csv').open('w',newline='') as f:
        wr=csv.writer(f); wr.writerow(['h_mm','max_relative_complex_error']); wr.writerows(rows)
    print(rows)
    assert all(x[1]/y[1]>3.8 for x,y in zip(rows,rows[1:]))
    rng=np.random.default_rng(977)
    m=b.data(20.,.025); y=m+.01*abs(m).max()*(rng.normal(size=(300,*m.shape))+1j*rng.normal(size=(300,*m.shape)))/np.sqrt(2)
    ref=fit(y,2.,b.W)
    for start in [(.08,-.8,.01),(.4,.9,.1)]:
        alt=fit(y,2.,b.W,start=start)
        print('initialization',start,'largest relative difference',abs(alt/ref-1).max())

if __name__=='__main__': main()
