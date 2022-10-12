from brainrender import Scene
from brainrender.actors.streamlines import Streamlines

import pandas as pd
import json
import glob
import re
from paper.figures import INSET # removes the scale bars?

###
# a few parameters / options
# a few parameters / options
interValue = True
saveScreen = True
animalName = "CR_transectomy3"
pltConnection = "UncutOnly" # "none" / "all" / "V1only"
streamExp = ["all"]  #["606100558"] list of experiments, or "all" if want them all
pltCut = "Cut" # "Cut" or "NoCut"
pltView = "side" # "coronal" or "side"
pltSlice = 0 # will slice the brain
dataPath = "C:/Users/Hamish/OneDrive - University College London/Paper/InjectionsA1_AllenBrainAtlas/"
screenFolder = dataPath + "screenshots"

###
# create scene
scene = Scene(atlas_name="allen_mouse_25um",inset=INSET,screenshots_folder=screenFolder)

# add brain regions
scene.add_brain_region("VISp","AUD", alpha=0.2)

# plot fibers
if pltConnection not in ["none"]:
    
    # if all get all regions
    if streamExp[0] in "all":
        expList = glob.glob(dataPath + "streamlines_*.json")
        streamExp = []
        for i in range(len(expList)):
            streamExp.append(re.search('streamlines_(.*).json', expList[i]).group(1))
    
    # do it one by one? There must be a faster and more elegant way
    for i in range(len(streamExp)):
        # get path
        if pltConnection == "all":
            streamlinePath = [dataPath + "streamlines_" + streamExp[i] + ".json"];
        else: 
            if pltConnection == "UncutOnly":
                streamlinePath = [dataPath + "processed/" + animalName + "_UnCutV1fibers_streamlines_" + streamExp[i] + ".json"]
            else:                 
                streamlinePath = [dataPath + "processed/" + animalName + "_UnCutV1fibers_streamlines_" + streamExp[i] + ".json",
                                   dataPath + "processed/" + animalName + "_cutV1fibers_streamlines_" + streamExp[i] + ".json"]
            
        # get color    
        if pltCut == "Cut":
            colors = ['darkred','salmon']
        else:
            colors = ['darkred','darkred']
        
        # plot streamline actors
        for i in range(len(streamlinePath)):
            with open(streamlinePath[i]) as json_data:
                     data = json.load(json_data)
            df = pd.json_normalize(data)
            if data["lines"] != []:
                # if not empty
                scene.add(Streamlines(df, color=colors[i], alpha=1))

# add cut
if pltCut == "Cut":
    cutPath = "D:/BrainSaw/" + animalName + "/downsampled_stacks/025_micron/brainreg_output/manual_segmentation/standard_space/regions/cut.obj";
    cut = scene.add(cutPath, color='black')
    scene.add_silhouette(cut, lw=2)

# add probe tracks?

# get camera right
if pltView == "coronal":
    cam = {
        "pos": (-36430, 0, -5700),
        "viewup": (0, -1, 0),
        "clippingRange": (30360, 64977),
        "focalPoint": (7319, 2861, -3942),
        "distance": 43901,
    }
elif pltView == "side":
    cam = {
        "pos": (11654, -32464, 21761),
        "viewup": (0, -1, -1),
        "clippingRange": (32024, 63229),
        "focalPoint": (7319, 2861, -3942),
        "distance": 43901,
    }
    
# render scene
if pltSlice:
    scene.slice("frontal")

scene.render(interactive=interValue, camera = cam, zoom = 1.5)
if saveScreen and not interValue:
    # deactivate interactive parameters in scene.render to be able to do that. Press "s" if doesn't work
    scene.screenshot(name="transectomy_" + animalName + "_" + streamExp[0] + "_Fibers-" + pltConnection + "_" + pltCut + "_" + pltView) 

    scene.close()