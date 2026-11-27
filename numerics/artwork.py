"""Typeset the measured layout and the two Monte Carlo maps as vector PDFs."""
import csv
from pathlib import Path
from reportlab.pdfgen.canvas import Canvas
from reportlab.lib.colors import HexColor,Color

BASE=Path(__file__).parent
def layout():
    c=Canvas(str(BASE/'layout.pdf'),pagesize=(450,92))
    c.setStrokeColor(HexColor('#334550'))
    c.line(28,48,421,48)
    for j,(kind,label) in enumerate([('heater','A'),('receiver','1'),('receiver','2'),('heater','B')]):
        x=55+112*j
        c.setFillColor(HexColor('#b85937' if kind=='heater' else '#286778'))
        if kind=='heater': c.rect(x-4,42,8,12,fill=1,stroke=0)
        else: c.circle(x,48,5,fill=1,stroke=0)
        c.setFillColorRGB(.1,.12,.14); c.setFont('Helvetica-Bold',10)
        c.drawCentredString(x,65,kind.capitalize()+' '+label)
        c.setFont('Helvetica',9);c.drawCentredString(x,27,['0','L','2L','3L'][j])
    c.setFont('Helvetica',9);c.drawCentredString(225,7,'Each heater is driven separately. The two receivers record each drive.')
    c.save()

def maps():
    with (BASE/'regimes.csv').open() as f: data=list(csv.DictReader(f))
    c=Canvas(str(BASE/'maps.pdf'),pagesize=(470,180))
    for p,(method,title) in enumerate([('ratio','Amplitude inversion'),('complex','Complex fit, gains unknown')]):
        x0=46+230*p;y0=31;cw=31;ch=22
        c.setFont('Helvetica-Bold',11);c.drawString(x0,159,title)
        for rec in data:
            x=sorted(set(float(r['Pe']) for r in data)).index(float(rec['Pe']))
            y=sorted(set(float(r['loss_number']) for r in data)).index(float(rec['loss_number']))
            value=float(rec[method])/200
            c.setFillColor(Color(.94-.83*value,.97-.57*value,.97-.50*value))
            c.rect(x0+x*cw,y0+y*ch,cw-1,ch-1,fill=1,stroke=0)
            c.setFont('Helvetica',9);c.setFillColorRGB(*((1,1,1) if value>.68 else (.08,.11,.12)))
            c.drawCentredString(x0+(x+.5)*cw,y0+(y+.5)*ch-3,f'{100*value:g}')
        c.setFillColorRGB(.1,.12,.14);c.setFont('Helvetica',9)
        for i,t in enumerate(['0.4','0.8','1.33','2','3']):c.drawCentredString(x0+(i+.5)*cw,19,t)
        for i,t in enumerate(['0.2','0.4','0.8','1.2','1.6']):c.drawRightString(x0-5,y0+(i+.5)*ch-3,t)
        c.drawCentredString(x0+2.5*cw,5,'Peclet number')
        c.saveState();c.translate(x0-28,y0+2.5*ch);c.rotate(90);c.drawCentredString(0,0,'Loss number');c.restoreState()
    c.save()
if __name__=='__main__': layout();maps()
