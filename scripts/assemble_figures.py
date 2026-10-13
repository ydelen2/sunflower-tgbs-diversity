#!/usr/bin/env python3
"""Assemble the multi-panel figures from the single-panel PNGs
written by the R steps."""

import argparse
import os
import sys

from PIL import Image, ImageDraw, ImageFont

LABEL_FONTS = [
    r"C:\Windows\Fonts\arialbd.ttf",
    "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
]

GAP = 60
MARGIN = 40
LABEL_PAD = 18
# Panel letter sizes in pixels, chosen so that the letters print at about
# 10 pt at the width each figure takes on the page.
LETTER_PX = {"structure": 120, "pca": 90, "fst": 95, "ld": 95, "qc": 100}


def load_font(size):
    for path in LABEL_FONTS:
        if os.path.exists(path):
            return ImageFont.truetype(path, size)
    try:
        import matplotlib

        path = os.path.join(
            os.path.dirname(matplotlib.__file__),
            "mpl-data", "fonts", "ttf", "DejaVuSans-Bold.ttf",
        )
        return ImageFont.truetype(path, size)
    except Exception:
        return ImageFont.load_default()


def scale_to_height(im, height):
    w = max(1, round(im.width * height / im.height))
    return im.resize((w, height), Image.LANCZOS)


def scale_to_width(im, width):
    h = max(1, round(im.height * width / im.width))
    return im.resize((width, h), Image.LANCZOS)


def label_band(font):
    ascent, descent = font.getmetrics()
    return ascent + descent + LABEL_PAD


def paste_panel(canvas, draw, im, x, y, letter, font):
    draw.text((x, y), letter, fill=(0, 0, 0), font=font)
    canvas.paste(im, (x, y + label_band(font)))


def assemble_figure1(src, out, row_height=1800):
    cv = scale_to_height(Image.open(os.path.join(src, "Fig3a_admixture_cv.png")), row_height)
    dapc = scale_to_height(Image.open(os.path.join(src, "Fig5_dapc.png")), row_height)
    font = load_font(LETTER_PX["structure"])
    band = label_band(font)

    top_width = cv.width + GAP + dapc.width
    bar = scale_to_width(Image.open(os.path.join(src, "Fig3c_admixture_K5.png")), top_width)

    width = top_width + 2 * MARGIN
    height = 2 * MARGIN + band + row_height + GAP + band + bar.height
    canvas = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(canvas)

    paste_panel(canvas, draw, cv, MARGIN, MARGIN, "a", font)
    paste_panel(canvas, draw, dapc, MARGIN + cv.width + GAP, MARGIN, "b", font)
    paste_panel(canvas, draw, bar, MARGIN, MARGIN + band + row_height + GAP, "c", font)

    canvas.save(out, dpi=(400, 400))
    return canvas.size


def assemble_figure2(src, out, row_height=2200):
    left = scale_to_height(Image.open(os.path.join(src, "Fig2b_pca_PC1PC2_status.png")), row_height)
    right = scale_to_height(Image.open(os.path.join(src, "Fig2e_pca_by_oil_content.png")), row_height)
    font = load_font(LETTER_PX["pca"])
    band = label_band(font)

    width = 2 * MARGIN + left.width + GAP + right.width
    height = 2 * MARGIN + band + row_height
    canvas = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(canvas)

    paste_panel(canvas, draw, left, MARGIN, MARGIN, "a", font)
    paste_panel(canvas, draw, right, MARGIN + left.width + GAP, MARGIN, "b", font)

    canvas.save(out, dpi=(400, 400))
    return canvas.size


def assemble_two_panels(src, out, left_name, right_name, row_height=2000, letter_px=95):
    """Two single-panel figures side by side, labelled a and b."""
    left = scale_to_height(Image.open(os.path.join(src, left_name)), row_height)
    right = scale_to_height(Image.open(os.path.join(src, right_name)), row_height)
    font = load_font(letter_px)
    band = label_band(font)
    width = 2 * MARGIN + left.width + GAP + right.width
    height = 2 * MARGIN + band + row_height
    canvas = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(canvas)
    paste_panel(canvas, draw, left, MARGIN, MARGIN, "a", font)
    paste_panel(canvas, draw, right, MARGIN + left.width + GAP, MARGIN, "b", font)
    canvas.save(out, dpi=(400, 400))
    return canvas.size


def assemble_figure3(src, out):
    """Pairwise differentiation (a) and genomic inbreeding (b) by cluster."""
    return assemble_two_panels(src, out, "Fig6_pairwise_fst.png", "Fig7_FROH_by_group.png",
                               letter_px=LETTER_PX["fst"])


def assemble_figure6(src, out):
    """LD decay on a linear axis (a) and on a logarithmic axis (b)."""
    return assemble_two_panels(src, out, "Fig9_ld_decay.png", "Fig9b_ld_decay_log.png",
                               letter_px=LETTER_PX["ld"])


def assemble_qc_panels(src, out, cell_height=1500):
    """Four quality-control histograms in a two by two grid, labelled a to d."""
    names = ["FigS2_sample_missingness.png", "FigS3_maf_spectrum.png",
             "FigS4_site_missing.png", "FigS5_site_het.png"]
    ims = [scale_to_height(Image.open(os.path.join(src, n)), cell_height) for n in names]
    cell_width = max(im.width for im in ims)
    font = load_font(LETTER_PX["qc"])
    band = label_band(font)
    width = 2 * MARGIN + 2 * cell_width + GAP
    height = 2 * MARGIN + 2 * (band + cell_height) + GAP
    canvas = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(canvas)
    for i, (im, letter) in enumerate(zip(ims, "abcd")):
        col, row = i % 2, i // 2
        x = MARGIN + col * (cell_width + GAP)
        y = MARGIN + row * (band + cell_height + GAP)
        paste_panel(canvas, draw, im, x, y, letter, font)
    canvas.save(out, dpi=(400, 400))
    return canvas.size


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    jobs = [
        ("Figure1_population_structure.png", assemble_figure1),
        ("Figure2_pca_status_and_oil.png", assemble_figure2),
        ("Figure3_fst_and_froh.png", assemble_figure3),
        ("Figure6_ld_decay_panels.png", assemble_figure6),
        ("FigureS4_quality_panels.png", assemble_qc_panels),
    ]
    for name, fn in jobs:
        dest = os.path.join(args.out, name)
        size = fn(args.src, dest)
        print(f"{name}\t{size[0]}x{size[1]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
