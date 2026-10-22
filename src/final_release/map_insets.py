"""Reuse the supplied PowerPoint locator without introducing result pixels."""
from pathlib import Path
import json
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image

def add_geographic_insets(fig,name,root):
    panels={'Fig2':[0,1,2],'Fig3':[0,1],'Fig4':[0],'Fig_S1':[0]}
    if name not in panels:return
    root=Path(root);out=root/'final_results';assets=out/'assets'
    im=Image.open(assets/'south_china_sea_reference.png')
    fig.canvas.draw()
    fig.savefig(assets/(name+'_without_insets.png'),dpi=250,bbox_inches='tight')
    fig.set_dpi(250);fig.canvas.draw()
    canvas=Image.fromarray(np.asarray(fig.canvas.buffer_rgba()).copy()).convert("RGB")
    boxes=[];files=[]
    for idx in panels[name]:
        ax=fig.axes[idx];bb=ax.get_window_extent();width=.20
        height=width*bb.width/bb.height*im.height/im.width
        if height>.32:width*=.32/height;height=.32
        # Crop this panel's own rendered southern region, retaining its exact colours.
        x0,y0=ax.transData.transform((105,18));x1,y1=ax.transData.transform((125,26))
        patch=canvas.crop((round(x0),round(canvas.height-y1),round(x1),round(canvas.height-y0)))
        cropfile=f'{name}_panel_{idx+1}_southern_crop.png';patch.save(assets/cropfile)
        composite=im.resize((744,992),Image.Resampling.LANCZOS)
        # Place the crop inside the retained reference frame, above the island symbols.
        bounds=(16,16,720,322)
        composite.paste(patch.resize((bounds[2]-bounds[0],bounds[3]-bounds[1]),Image.Resampling.LANCZOS),bounds[:2])
        imagefile=f'{name}_panel_{idx+1}_inset.png';composite.save(assets/imagefile);files.append(imagefile)
        ia=ax.inset_axes([1-width-.025,.025,width,height],zorder=20)
        ia.imshow(composite,interpolation='lanczos');ia.set_axis_off()
        boxes.append(ia)
    fig.canvas.draw();renderer=fig.canvas.get_renderer()
    crop=fig.get_tightbbox(renderer).padded(.1)
    positions=[]
    for ia,imagefile in zip(boxes,files):
        bb=ia.get_window_extent(renderer).transformed(fig.dpi_scale_trans.inverted())
        positions.append({'image':imagefile,'left':(bb.x0-crop.x0)/crop.width,'top':(crop.y1-bb.y1)/crop.height,'width':bb.width/crop.width,'height':bb.height/crop.height})
    (assets/(name+'_inset_positions.json')).write_text(json.dumps(positions,indent=2))
