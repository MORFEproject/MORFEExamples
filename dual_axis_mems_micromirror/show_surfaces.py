"""ParaView view of the numbered boundary surfaces written by the notebook.

Colours every boundary surface of the device layer and writes its `surface_id`
on it. Run it from this folder, after section 1 of the notebook, with

    paraview --script=show_surfaces.py

(on macOS: /Applications/ParaView-<version>.app/Contents/MacOS/paraview), or from
ParaView's Python Shell with Run Script. With `pvpython show_surfaces.py out.png`
it saves a top view instead.

The numbers are ParaView selection labels: clicking in the view clears them, and
re-running the script brings them back. Hover Cells On (toolbar) also shows the
`surface_id` of the face under the cursor. The surfaces currently clamped in the
notebook carry `clamped = 1`.
"""
import os
import sys

from paraview.simple import (ColorBy, GetActiveViewOrCreate, GetColorTransferFunction,
                             QuerySelect, Render, ResetCamera, SaveScreenshot,
                             SetActiveSource, Show, XMLUnstructuredGridReader)

try:
    here = os.path.dirname(os.path.abspath(__file__))
except NameError:  # ParaView's Run Script does not always set __file__
    here = os.getcwd()
pv_dir = os.path.join(here, "results", "paraview")

surfaces = XMLUnstructuredGridReader(registrationName="surfaces",
                                     FileName=[os.path.join(pv_dir, "surfaces.vtu")])
labels = XMLUnstructuredGridReader(registrationName="surface_labels",
                                   FileName=[os.path.join(pv_dir, "surface_labels.vtu")])

view = GetActiveViewOrCreate("RenderView")

# `color_key` shuffles the ids so that neighbouring surfaces get distinct colours.
s = Show(surfaces, view)
s.SetRepresentationType("Surface With Edges")
s.EdgeColor = [0.25, 0.25, 0.25]
ColorBy(s, ("CELLS", "color_key"))
GetColorTransferFunction("color_key").ApplyPreset("Rainbow Uniform", True)
s.SetScalarBarVisibility(view, False)

l = Show(labels, view)
l.SetRepresentationType("Points")
l.PointSize = 1
l.SelectionPointFieldDataArrayName = "surface_id"
l.SelectionPointLabelVisibility = 1
l.SelectionPointLabelFontSize = 18
l.SelectionPointLabelBold = 1
l.SelectionPointLabelColor = [0.0, 0.0, 0.0]
l.SelectionColor = [0.0, 0.0, 0.0]
SetActiveSource(labels)
QuerySelect(QueryString="surface_id >= 0", FieldType="POINT", InsideOut=0)

# Look down onto the handle side (z = 45), where the bonded face is.
view.CameraPosition = [0.0, 0.0, 6000.0]
view.CameraFocalPoint = [0.0, 0.0, 0.0]
view.CameraViewUp = [0.0, 1.0, 0.0]
ResetCamera(view)
Render(view)

if len(sys.argv) > 1:
    view.UseColorPaletteForBackground = 0
    view.Background = [1.0, 1.0, 1.0]
    view.OrientationAxesVisibility = 1
    view.ViewSize = [1800, 1400]
    SaveScreenshot(sys.argv[1], view, ImageResolution=[1800, 1400])
