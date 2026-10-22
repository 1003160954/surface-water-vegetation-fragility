"""Reuse the supplied PowerPoint locator without introducing result pixels."""
from pathlib import Path
import json
import matplotlib.pyplot as plt
from PIL import Image

def add_geographic_insets(fig,name,root):
    panels={'Fig2':[0,1,2],'Fig3':[0,1],'Fig4':[0],'Fig_S1':[0]}
    if name not in panels:return
    root=Path(root);out=root/'final_results';assets=out/'assets'
    im=Image.open(assets/'south_china_sea_reference.png')
    fig.canvas.draw()
    fig.savefig(assets/(name+'_without_insets.png'),dpi=250,bbox_inches='tight')
    boxes=[]
    for idx in panels[name]:
        ax=fig.axes[idx];bb=ax.get_window_extent();width=.20
        height=width*bb.width/bb.height*im.height/im.width
        if height>.32:width*=.32/height;height=.32
        ia=ax.inset_axes([1-width-.025,.025,width,height],zorder=20)
        ia.imshow(im,interpolation='lanczos');ia.set_axis_off()
        boxes.append(ia)
    fig.canvas.draw();renderer=fig.canvas.get_renderer()
    crop=fig.get_tightbbox(renderer).padded(.1)
    positions=[]
    for ia in boxes:
        bb=ia.get_window_extent(renderer).transformed(fig.dpi_scale_trans.inverted())
        positions.append({'left':(bb.x0-crop.x0)/crop.width,'top':(crop.y1-bb.y1)/crop.height,'width':bb.width/crop.width,'height':bb.height/crop.height})
    (assets/(name+'_inset_positions.json')).write_text(json.dumps(positions,indent=2))
