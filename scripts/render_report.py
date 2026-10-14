"""Optional report renderer: pip install reportlab; python scripts/render_report.py."""
from pathlib import Path
import json
import csv
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
import matplotlib
fontroot=Path(matplotlib.get_data_path())/"fonts/ttf"
pdfmetrics.registerFont(TTFont("DejaVu",str(fontroot/"DejaVuSans.ttf")))
pdfmetrics.registerFont(TTFont("DejaVu-Bold",str(fontroot/"DejaVuSans-Bold.ttf")))
pdfmetrics.registerFontFamily("DejaVu",normal="DejaVu",bold="DejaVu-Bold",italic="DejaVu",boldItalic="DejaVu-Bold")
from reportlab.lib.colors import HexColor, white
from reportlab.platypus import Paragraph, Table, TableStyle
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.utils import ImageReader

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs'
m=json.loads((OUT/'metadata.json').read_text())
q=m['qc']
W,H=842,595
NAVY=HexColor('#153449'); TEAL=HexColor('#247c8c'); GRAY=HexColor('#4b626e')
path=OUT/'QUEST_Gravity_Technical_Brief.pdf'
c=canvas.Canvas(str(path),pagesize=(W,H))
c.setTitle('QUEST airborne gravity - Reproducible 2D analysis')
c.setAuthor('AI-assisted portfolio project; source: Geoscience BC / Sander Geophysics')
styles={
 'body':ParagraphStyle('body',fontName='DejaVu',fontSize=10,leading=14,textColor=NAVY),
 'small':ParagraphStyle('small',fontName='DejaVu',fontSize=8.4,leading=11,textColor=GRAY),
 'head':ParagraphStyle('head',fontName='DejaVu-Bold',fontSize=12,leading=16,textColor=NAVY),
}
def para(text,x,y,w=750,style='body'):
 p=Paragraph(text,styles[style]);_,h=p.wrap(w,H);p.drawOn(c,x,y-h);return y-h

def base(number,title,subtitle):
 c.setFillColor(NAVY);c.rect(0,H-9,W,9,fill=1,stroke=0)
 c.setFont('DejaVu-Bold',23);c.drawString(40,H-47,title)
 para(subtitle,40,H-62,760,'small')
 c.setStrokeColor(HexColor('#dce5e9'));c.line(40,35,W-40,35)
 c.setFillColor(GRAY);c.setFont('DejaVu',8)
 c.drawString(40,21,'QUEST / Geoscience BC Report 2008-08 / Sander Geophysics  |  Reproducible Python case study')
 c.drawRightString(W-40,21,f'{number} / 5')

def fig(name,x,y,w):
 p=OUT/'figures'/name
 im=ImageReader(str(p));iw,ih=im.getSize();h=w*ih/iw
 c.drawImage(im,x,y-h,width=w,height=h,mask='auto');return y-h

def table(rows,x,y,widths):
 formatted=[[Paragraph(str(v),styles['small']) for v in row] for row in rows]
 t=Table(formatted,colWidths=widths,hAlign='LEFT')
 t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),HexColor('#e4eff2')),
                       ('ROWBACKGROUNDS',(0,1),(-1,-1),[white,HexColor('#f4f7f9')]),
                       ('VALIGN',(0,0),(-1,-1),'TOP'),('BOTTOMPADDING',(0,0),(-1,-1),8),
                       ('TOPPADDING',(0,0),(-1,-1),8),('LEFTPADDING',(0,0),(-1,-1),8),
                       ('RIGHTPADDING',(0,0),(-1,-1),8)]))
 _,h=t.wrap(sum(widths),H);t.drawOn(c,x,y-h);return y-h

base(1,'QUEST airborne gravity','2D anomaly mapping, metadata auditing and held-line interpolation validation | Original QUEST survey, British Columbia')
y=table([['OBSERVATIONS','LINE SECTIONS','LINEAR CV RMSE','SUPPORTED TEST BLOCKS'],
         [f"{q['rows_parsed']:,}",str(q['line_sections']),'1.912 mGal','98.8%']],40,500,[190]*4)
y=fig('01_gravity_atlas.png',70,y-8,700)
para('<b>Outcome.</b> The supplier fields were converted into GIS-ready rasters, and a separate line-data interpolation was tested with whole planned lines withheld. The 250 m source-grid spacing is not a claim of 250 m geophysical resolution.',40,92,760)
c.showPage()

base(2,'Audit before interpretation','The release contains conflicting metadata. Each decision is recorded and checked; the original source files remain unchanged.')
rows=[['Issue','Decision and measured evidence'],
 ['CRS conflict','README/XYZ specify NAD83 / UTM 10N; report table says zone 32N. Use EPSG:26910. A 10,496-point coordinate check gives a median discrepancy of 1.303 m, supporting the zone selection. This is not positional accuracy.'],
 ['File format','Released grids are binary Geosoft files despite an ASCII GXF label in the report. A guarded reader accepts only the observed encoding. Actual grid-versus-line sampling checks values and orientation.'],
 ['Gravity channel labels','Some final channel descriptions incorrectly suggest isostatic correction on every field. Grid correspondence and the isostatic-channel identity support treating the three products separately. Do not apply corrections twice.'],
 ['Correction history','The line table says simple Bouguer, while the narrative includes terrain and curvature. GRVBGL100 is labelled as the supplier Bouguer channel without claiming a newly computed correction decomposition.'],
 ['Time and missing values','43,663 records have SECOND >= 86,400 and require rollover interpretation for time-based work. No missing numeric values or duplicate flight/line/time records were found. Spatial analysis retains these observations.']]
y=table(rows,40,500,[140,620])
y=para('<b>Source consistency checks</b>',40,y-22,760,'head')
cons=list(csv.DictReader((OUT/'tables/supplier_grid_consistency.csv').open()))
rows=[['Supplier grid / line channel','RMSE of related products (mGal)']]+[[r['file']+' / '+r['line_channel'],f"{float(r['rmse_mgal']):.5f}"] for r in cons]
y=table(rows,40,y-10,[500,260])
para('These tiny discrepancies verify correspondence of related products. They are not independent validation against ground truth. Full evidence and unresolved details are in docs/METADATA_ISSUES.md.',40,y-14,760,'small')
c.showPage()

base(3,'Line observations to supported rasters','A second workflow uses GRVBGL100, preserving a clear distinction from the supplier\'s final spatially filtered grids.')
y=fig('03_line_processing.png',40,505,760)
y=para('<b>Method.</b> Average samples in 500 m XY cells separately for each exact line section. Use traverse/infill lines for linear triangulation; exclude control lines. Export a 1 km node raster. Mask outside the convex hull, beyond 3 km from data or inside triangles with edges longer than 7 km.',40,y-10,760)
para('<b>Interpretation.</b> Distance to data is a support diagnostic, not a formal uncertainty estimate. These masks limit gap bridging. Linear interpolation can leave triangular features; a visually detailed surface is not evidence of finer geophysical resolving power.',40,y-12,760)
c.showPage()

base(4,'Whole-line validation','Five interleaved folds; all sections and reflights of a planned line remain together. Control lines are excluded from training and testing.')
y=fig('04_validation.png',40,500,760)
rows=[['Interpolator','RMSE (mGal)','MAE (mGal)','Bias (mGal)','Coverage']]
for r in csv.DictReader((OUT/'tables/cv_summary.csv').open()):
 rows.append([r['method'],f"{float(r['rmse_mgal']):.3f}",f"{float(r['mae_mgal']):.3f}",f"{float(r['bias_mgal']):.3f}",f"{float(r['coverage_fraction']):.1%}"])
y=table(rows,40,y-5,[240,130,130,130,130])
para('<b>Limit.</b> The errors describe interpolation of already corrected and levelled data. Supplier-wide levelling, nearby infill lines and spatial correlation remain. No independent instrument accuracy or operational navigation performance is inferred. Linear was fixed in advance; the small difference versus IDW is not a statistical-significance claim.',40,y-14,760,'small')
c.showPage()

base(5,'Spatial analysis and the next experiment','Exploratory regional/residual separation and horizontal gradients connect this study to questions about map structure.')
y=fig('02_spatial_analysis.png',80,505,680)
y=para('<b>Scale and edges.</b> The residual subtracts a Gaussian regional component with sigma = 20 km and requires 95% kernel support. On a common 11,969-cell area, residual standard deviation changes from 6.20 to 7.99 to 9.22 mGal for sigma = 10, 15 and 20 km. No cells support sigma = 30 km at the chosen threshold.',40,y-6,760)
y=para('<b>Navigation relevance.</b> Gradients describe spatial variation, not unique position observability. Before a navigation simulation, agree on the gravity functional, observation height, terrain/long-wavelength treatment, independent validation data and realistic sensor/map errors.',40,y-10,760)
para('<b>Reproduce:</b> README.md, config.json, requirements-lock.txt and outputs/metadata.json. Eight numerical/GIS tests passed. AI-assisted initial implementation; candidate review and independent rerun are recommended before presentation. Source: geosciencebc.com/major-projects/quest/.',40,y-10,760,'small')
c.save()
print(path)
