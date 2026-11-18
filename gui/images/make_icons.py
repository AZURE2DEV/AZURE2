#!/usr/bin/env python3
"""Render the GUI icon files from the SVG logos (2026-11-06).

    python3 gui/images/make_icons.py        # needs rsvg-convert and Pillow

azure2_logo.svg            -> azure2_icon.png (256 px; the default window icon, logo button, About box,
                              Linux AppImage), azure2.ico (16-256 px), azure2.icns (16-1024 px)
azure2_logo_polarized.svg  -> azure2_icon_polarized.png (256 px), which the GUI shows instead while the
                              project holds a polarization segment (AZURESetup::updateLogo)
"""
import io, os, subprocess
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))


def render(svg, px):
    png = subprocess.run(["rsvg-convert", "-w", str(px), "-h", str(px), os.path.join(HERE, svg)],
                         check=True, capture_output=True).stdout
    return Image.open(io.BytesIO(png)).convert("RGBA")


render("azure2_logo.svg", 256).save(os.path.join(HERE, "azure2_icon.png"), optimize=True)
render("azure2_logo_polarized.svg", 256).save(os.path.join(HERE, "azure2_icon_polarized.png"), optimize=True)

ico_sizes = [16, 22, 24, 32, 40, 48, 64, 96, 128, 256]
ico = {s: render("azure2_logo.svg", s) for s in ico_sizes}
ico[256].save(os.path.join(HERE, "azure2.ico"), format="ICO", sizes=[(s, s) for s in ico_sizes],
              append_images=[ico[s] for s in ico_sizes if s != 256])

icns_sizes = [16, 32, 64, 128, 256, 512, 1024]
icns = {s: render("azure2_logo.svg", s) for s in icns_sizes}
icns[1024].save(os.path.join(HERE, "azure2.icns"), format="ICNS",
                append_images=[icns[s] for s in icns_sizes if s != 1024])
print("wrote azure2_icon.png, azure2_icon_polarized.png, azure2.ico, azure2.icns")
