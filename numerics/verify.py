"""Harmonic heat transport; mm and seconds. Run with Python 3 + NumPy."""
from pathlib import Path
import csv
import numpy as np

D, V, LOSS, L = .15, .10, .03, 2.
W = np.array([.025,.12])

def forward(w,v,margin,step):
    """Zero-temperature ends, central finite differences, point heater loads."""
    n=round((2*margin+3*L)/step)-1
    diag=np.full(n,2*D/step**2+LOSS+1j*w)
    lower=-D/step**2-v/(2*step)
    upper=-D/step**2+v/(2*step)
    rhs=np.zeros((n,2),complex)
    def index(x): return round((x+margin)/step)-1
    rhs[index(0.),0]=rhs[index(3*L),1]=1/step
    for i in range(1,n):
        mult=lower/diag[i-1]
        diag[i]-=mult*upper
        rhs[i]-=mult*rhs[i-1]
    rhs[-1]/=diag[-1]
    for i in reversed(range(n-1)):
        rhs[i]=(rhs[i]-upper*rhs[i+1])/diag[i]
    return rhs[[index(L),index(2*L)]]

def infer(m):
    # Leading axes: flow direction, angular frequency, receiver, heater.
    y=np.log(abs(m))
    a=(y[:,:,0,0]+y[:,:,1,1]-y[:,:,1,0]-y[:,:,0,1]).mean(axis=0)/(2*L)
    d=.5*np.sqrt(np.diff((W/a)**2)[0]/np.diff(a*a)[0])
    beta=(y[0,:,1,:]-y[0,:,0,:]-y[1,:,1,:]+y[1,:,0,:]).mean()/L
    loss=d*(a[0]**2-(W[0]/(2*d*a[0]))**2-beta**2/4)
    return np.array([d,d*beta,loss])

def data(margin,step):
    raw=np.array([[forward(w,v,margin,step) for w in W] for v in [V,-V]])
    # Different complex receiver gains and lags; heater amplitudes change on reversal.
    sensors=np.array([1.3,.72])[None,:]/(1+1j*W[:,None]*np.array([3.,9.]))
    heaters=np.array([[.7,1.4],[1.1,.8]])[:,None,None,:]
    return raw*sensors[None,:,:,None]*heaters

def main():
    truth=np.array([D,V,LOSS])
    rng=np.random.default_rng(51425)
    rows=[]
    for margin,step in [(2.,.05),(6.,.05),(20.,.10),(20.,.05),(20.,.025)]:
        m=data(margin,step)
        rows.append([margin,step,'clean',*infer(m)])
    np.testing.assert_allclose(infer(m),truth,rtol=.001)
    for kind in ['relative','absolute']:
        vals=[]
        for _ in range(3000):
            scale=abs(m) if kind=='relative' else abs(m).max()
            perturb=.01*scale*(rng.normal(size=m.shape)+1j*rng.normal(size=m.shape))/np.sqrt(2)
            with np.errstate(invalid='ignore'):
                vals.append(infer(m+perturb))
        vals=np.array(vals)
        valid=np.isfinite(vals).all(axis=1)
        print(kind,'finite',valid.sum(),'negative loss',(vals[:,2]<0).sum())
        for q,entry in zip([5,50,95],np.percentile(vals[valid],[5,50,95],axis=0)):
            rows.append([20.,.025,kind+'_p'+str(q),*entry])
    with Path(__file__).with_name('results.csv').open('w',newline='') as f:
        wr=csv.writer(f);wr.writerow(['margin_mm','step_mm','scenario','D_mm2_s','v_mm_s','loss_s_inv']);wr.writerows(rows)
    for row in rows: print(*row)

if __name__=='__main__': main()
