"""Paired simulations. Python >=3.10, NumPy >=2.0. Time in s, length in mm."""
from pathlib import Path
import csv
import numpy as np
import verify as base

OUT=Path(__file__).resolve().parent
truth=np.array([.15,.1,.03])
direction=np.array([1.,-1.])[None,:,None,None]
side=np.array([-1.,1.])[None,None,None,:]
def residual(p,y,L,w):
    diffus=np.exp(p[:,0])
    beta=p[:,1]
    loss=np.exp(p[:,2])
    k=np.sqrt(beta[:,None]**2/4+(loss[:,None]+1j*w)/diffus[:,None])
    gains=p[:,3:5]+1j*p[:,5:7]
    logq=gains[:,None,:,None]+direction*beta[:,None,None,None]*L/2+side*k[:,None,:,None]*L
    ratio=np.exp(logq.real.clip(-55,55)+1j*logq.imag)
    err=(y[:,:,:,1,:]-ratio*y[:,:,:,0,:])/np.sqrt(1+abs(ratio)**2)
    err=err.reshape(len(p),-1)
    return np.column_stack((err.real,err.imag))

def fit(y,L,w,calibrated=False,start=(.22,0.,.045)):
    n=len(y)
    p=np.tile([np.log(start[0]),start[1],np.log(start[2]),0.,0.,0.,0.],(n,1))
    if calibrated:
        ratio=(.72/1.3)*(1+3j*w)/(1+9j*w)
        p[:,3:5]=np.log(abs(ratio)); p[:,5:7]=np.angle(ratio)
    active=3 if calibrated else 7
    lam=np.full(n,.01)
    for _ in range(55):
        res=residual(p,y,L,w)
        eps=1.e-5
        jac=np.stack([(residual(p+dp,y,L,w)-res)/eps for dp in eps*np.eye(7)[:active]],axis=2)
        jt=jac.swapaxes(1,2)
        h=jt@jac
        norm=np.maximum(h.diagonal(axis1=1,axis2=2),1.e-9)
        step=np.linalg.solve(h+np.eye(active)*norm[:,None,:]*lam[:,None,None],jt@res[:,:,None])[:,:,0]
        trial=p.copy()
        trial[:,:active]-=np.clip(step,-1,1)
        trial[:,[0,2]]=trial[:,[0,2]].clip(-12,3)
        ok=(residual(trial,y,L,w)**2).sum(1)<(res**2).sum(1)
        p=np.where(ok[:,None],trial,p)
        lam=np.where(ok,np.maximum(1.e-8,lam/3),np.minimum(1.e8,lam*5))
    d=np.exp(p[:,0])
    return np.column_stack((d,d*p[:,1],np.exp(p[:,2])))

def experiment(L,w,noise,N,seed,*,reverse=1.,drift=1.,mode=0.,velocity=.1,loss0=.03):
    base.L,base.W=L,w
    def solve(loss):
        base.LOSS=loss
        try:
            return np.array([[base.forward(om,v,20.,.025) for om in w]
                             for v in [velocity,-velocity*reverse]])
        finally:
            base.LOSS=.03
    response=solve(loss0)
    if mode: response+=mode*solve(loss0+.15*(np.pi/3)**2)
    g=np.array([1.3,.72])[None,:]/(1+1j*w[:,None]*np.array([3.,9.]))
    m=response*g[None,:,:,None]*np.array([[.7,1.4],[1.1,.8]])[:,None,None,:]
    m[0,:,1]*=drift
    rng=np.random.default_rng(seed)
    y=m+noise*(rng.normal(size=(N,*m.shape))+1j*rng.normal(size=(N,*m.shape)))/np.sqrt(2)
    with np.errstate(invalid='ignore',divide='ignore'):
        closed=np.array([base.infer(obs) for obs in y])
    return m,[closed,fit(y,L,w),fit(y,L,w,True)]

def main():
    w=np.array([.025,.12])
    base.L,base.W=2.,w
    ref=abs(base.data(20.,.025)).max()
    out=[]
    for case,L,ww,frac,rev,drift,mode in [
        ('base',2,w,.01,1,1,0),('quiet',2,w,.002,1,1,0),
        ('short',1,w,.01,1,1,0),('medium',1.5,w,.01,1,1,0),
        ('long',3,w,.01,1,1,0),('lowpair',2,np.array([.025,.06]),.01,1,1,0),
        ('highpair',2,np.array([.025,.20]),.01,1,1,0),
        ('unequal_reverse',2,w,.002,.9,1,0),('gain_drift',2,w,.002,1,1.05,0),
        ('transverse',2,w,.002,1,1,.35)]:
        m,results=experiment(L,ww,ref*frac,300,977,reverse=rev,drift=drift,mode=mode)
        for method,vals in zip(['ratio','complex','calibrated'],results):
            physical=np.isfinite(vals).all(1)&(vals[:,0]>0)&(vals[:,2]>=0)
            success=physical&(abs(vals/truth-1)<.1).all(1)
            med=np.median(vals[physical],axis=0)
            out.append([case,L,*ww,ref*frac,abs(m).min()/(ref*frac),method,
                        physical.sum(),success.sum(),*med,*np.median(abs(vals[physical]/truth-1),axis=0)])
    with (OUT/'benchmark.csv').open('w',newline='') as f:
        wr=csv.writer(f);wr.writerow('case L w1 w2 noise minSNR method physical within10 medD medV medLoss MdARE_D MdARE_v MdARE_loss'.split());wr.writerows(out)
    print('\n'.join(f'{r[0]:18} {r[6]:7} {r[5]:5.2f} {r[7]:3} {r[8]:3} '+str(np.round(r[9:12],4)) for r in out))
if __name__=='__main__':main()



