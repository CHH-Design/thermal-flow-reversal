"""Equal absolute noise floor throughout, 200 paired trials per cell."""
from pathlib import Path
import csv
import numpy as np
import verify as b
from study import experiment

def main():
    w=np.array([.025,.12])
    b.L,b.W=2.,w
    noise=.01*abs(b.data(20.,.025)).max()
    with Path(__file__).with_name('regimes.csv').open('w',newline='') as f:
        out=csv.writer(f); out.writerow(['Pe','loss_number','minSNR','ratio','complex','calibrated'])
        for pe in [.4,.8,1.33,2.,3.]:
            for loss_number in [.2,.4,.8,1.2,1.6]:
                v=pe*.15/2; loss=loss_number*.15/4
                m,est=experiment(2.,w,noise,200,1993,velocity=v,loss0=loss)
                counts=[(abs(z/np.array([.15,v,loss])-1)<.1).all(1).sum() for z in est]
                out.writerow([pe,loss_number,abs(m).min()/noise,*counts])

if __name__=='__main__': main()
