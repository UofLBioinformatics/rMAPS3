import os
import re
import shutil
import csv
import atexit

from pyx import *

pdfium = None
Image = None


class RmapsUnicodeEngine(text.UnicodeEngine):
    """Unicode text engine that tolerates legacy PyX text attributes."""

    def text_pt(self,
                x_pt,
                y_pt,
                text_value,
                textattrs=[],
                texmessages=[],
                fontmap=None,
                singlecharmode=False):
        merged = attr.mergeattrs(textattrs)
        supported = attr.getattrs(merged, [trafo.trafo_pt, style.fillstyle])
        return super().text_pt(
            x_pt,
            y_pt,
            text_value,
            supported,
            texmessages=texmessages,
            fontmap=fontmap,
            singlecharmode=singlecharmode,
        )


def use_unicode_text_engine(logger=None):
    """Use a non-TeX PyX text engine for environments where TeX is blocked."""
    try:
        text.set(RmapsUnicodeEngine)
    except Exception as exc:
        if logger:
            logger.warning("PyX Unicode text engine unavailable: %s", exc)
        return False
    return True


def boxes(xS, width, scale, boxY, box_height, splice_offset):
    xE = xS + width
    rect = path.rect(xS * scale, boxY * scale, width * scale,
                     box_height * scale)
    r_line = [
        path.line(xS * scale, (boxY + box_height / 4.0) * scale, xE * scale,
                  (boxY + box_height / 4.0) * scale),
        path.line(xS * scale, (boxY + box_height / 2.0) * scale, xE * scale,
                  (boxY + box_height / 2.0) * scale),
        path.line(xS * scale, (boxY + 3 * box_height / 4.0) * scale,
                  xE * scale, (boxY + 3 * box_height / 4.0) * scale)
    ]
    r_ss = path.line((xS + splice_offset) * scale, boxY * scale,
                     (xS + splice_offset) * scale, (boxY + box_height) * scale)
    return rect, r_line, r_ss


def title_and_legend(c, scale, xE, divider_gap, boxY, box_height, exon_width,
                     intron_width, map_name, nu, nd, nb):
    mls = style.linewidth.THIck
    bgmls = style.linewidth.Thick
    upColor = color.rgb.red
    dnColor = color.rgb.blue
    bgColor = color.rgb.black
    largeLab = text.size.huge
    titleLab = text.size.Huge
    ldash = style.linestyle.dashed
    #c.text((xE+gap)*scale, (boxY+box_height+45)*scale, motifLabel + " Motif MAP", [titleLab, text.halign.boxcenter]);
    c.text((xE + divider_gap) * scale, (boxY + box_height + 85) * scale,
           "Motif MAP: " + map_name.replace('_', ''),
           [titleLab, text.halign.boxcenter])

    c.stroke(
        path.line(
            (xE + divider_gap * 2) * scale, (boxY + box_height + 25) * scale,
            (xE + divider_gap * 3) * scale, (boxY + box_height + 25) * scale),
        [upColor, mls])
    #c.text((xE+gap*4)*scale, (boxY+box_height+20)*scale, "Upregulated ("+str(totalExonCount['up'])+')', [upColor,largeLab]);
    c.text((xE + divider_gap * 4) * scale, (boxY + box_height + 20) * scale,
           "Upregulated(" + str(nu) + ")", [upColor, largeLab])
    c.stroke(
        path.line((xE + divider_gap + divider_gap + intron_width * 2 / 3 +
                   5 * divider_gap) * scale, (boxY + box_height + 25) * scale,
                  (xE + divider_gap + divider_gap + intron_width * 2 / 3 +
                   divider_gap + 5 * divider_gap) * scale,
                  (boxY + box_height + 25) * scale), [dnColor, mls])
    #c.text((xE+gap+gap+iW*2/3+gap+gap)*scale, (boxY+box_height+20)*scale, "Downregulated ("+str(totalExonCount['dn'])+')', [dnColor,largeLab]);
    c.text((xE + divider_gap + divider_gap + intron_width * 2 / 3 +
            divider_gap + 5 * divider_gap) * scale,
           (boxY + box_height + 20) * scale, "Downregulated(" + str(nd) + ")",
           [dnColor, largeLab])
    c.stroke(
        path.line((xE + divider_gap + intron_width + divider_gap * 18) * scale,
                  (boxY + box_height + 25) * scale,
                  (xE + divider_gap + intron_width + divider_gap * 19) * scale,
                  (boxY + box_height + 25) * scale), [bgColor, bgmls])
    #c.text((xE+gap+iW+gap*11)*scale, (boxY+box_height+20)*scale, "Background ("+str(totalExonCount['bg'])+')', [bgColor,largeLab]);
    c.text((xE + divider_gap + intron_width + divider_gap * 20) * scale,
           (boxY + box_height + 20) * scale, "Background(" + str(nb) + ")",
           [bgColor, largeLab])

    c.stroke(
        path.line(
            (xE + divider_gap * 2) * scale, (boxY + box_height + 55) * scale,
            (xE + divider_gap * 3) * scale, (boxY + box_height + 55) * scale),
        [upColor, ldash])
    c.text((xE + divider_gap * 4) * scale, (boxY + box_height + 50) * scale,
           "-log(pVal) up vs. bg", [upColor, largeLab])
    c.stroke(
        path.line((xE + 8 * divider_gap + intron_width * 2 / 3) * scale,
                  (boxY + box_height + 55) * scale,
                  (xE + 8 * divider_gap + intron_width * 2 / 3 + divider_gap) *
                  scale, (boxY + box_height + 55) * scale), [dnColor, ldash])
    c.text((xE + 8 * divider_gap + intron_width * 2 / 3 + divider_gap +
            divider_gap) * scale, (boxY + box_height + 50) * scale,
           "-log(pVal) dn vs. bg", [dnColor, largeLab])


def make_label(max):
    return [
        '0.0',
        str("%.3f" % (max * 0.25)),
        str("%.3f" % (max * 0.5)),
        str("%.3f" % (max * 0.75)),
        str("%.3f" % float(max))
    ]


def draw_y_axis(c, scale, yLab, xS, boxY, box_height):
    yl_offset = 5
    largeLab = text.size.huge
    yLabAtt = [largeLab, text.halign.boxright, text.valign.middle]
    # y-axis
    nYDivisions = float(len(yLab) - 1)
    for i in range(len(yLab)):
        c.text((xS - yl_offset) * scale,
               (boxY + box_height * i / nYDivisions) * scale, yLab[i], yLabAtt)
    ## y-label
    c.text((xS - yl_offset * 15) * scale, (boxY + box_height / 2.0) * scale,
           "Motif Score (mean)",
           [largeLab, text.halign.boxcenter,
            trafo.rotate(90)])


def draw_yp_axis(c, scale, ypLab, yp_y_offset, xE, boxY, box_height,
                 divider_gap):
    yl_offset = 5
    largeLab = text.size.huge
    ypLabAtt = [largeLab, text.halign.boxleft, text.valign.middle]
    # yp-axis
    nYPDivisions = float(len(ypLab) - 1)
    for i in range(len(ypLab)):
        c.text((xE - 5 * yl_offset - divider_gap) * scale,
               (boxY + box_height * i / nYPDivisions + yp_y_offset) * scale,
               ypLab[i], ypLabAtt)
    ## y-label
    c.text((xE + yl_offset * 6) * scale,
           (boxY + box_height / 2.0 + yp_y_offset) * scale,
           "Negative log10(pValue)",
           [largeLab, text.halign.boxcenter,
            trafo.rotate(270)])


def draw_x_axis_segment(c, scale, xS, boxY, labels, units_per_bp):
    xl_offset = 23
    largeLab = text.size.huge
    for i in range(len(labels)):
        label = labels[i]
        t = str(label)
        offset = (label - labels[0]) * units_per_bp
        styles = [largeLab] if i < len(labels) - 1 else [
            largeLab, text.halign.boxright
        ]
        c.text((xS + offset) * scale, (boxY - xl_offset) * scale, t, styles)


def export_canvas_outputs(canv,
                          pdf_path,
                          png_path,
                          png_resolution=100,
                          logger=None):
    """
    Export a PyX canvas with robust fallbacks:
    1) PDF via PyX
    2) PNG via Ghostscript-backed writeGSfile
    3) Python fallback via pypdfium2 (render first PDF page to PNG)
    """
    pdf_ok = False
    png_ok = False

    global pdfium, Image

    try:
        canv.writePDFfile(pdf_path)
        pdf_ok = True
    except Exception as exc:
        if logger:
            logger.debug("Native PDF export failed for %s: %s", pdf_path, exc)

    if pdf_ok and pdfium is None:
        try:
            import pypdfium2 as pdfium_module
            pdfium = pdfium_module
        except Exception:
            pdfium = None

    if pdf_ok and pdfium is not None:
        try:
            doc = pdfium.PdfDocument(pdf_path)
            page = doc[0]
            scale = max(1.0, float(png_resolution) / 72.0)
            bitmap = page.render(scale=scale)
            image = bitmap.to_pil()
            image.save(png_path, format="PNG")
            png_ok = True
            doc.close()
        except Exception as exc:
            if logger:
                logger.debug("PDFium PNG fallback failed for %s: %s",
                             png_path, exc)

    if not png_ok:
        try:
            tmpdir = os.path.dirname(pdf_path) or "."
            gs_tmp = os.path.join(tmpdir, "jj.png")
            canv.writeGSfile(gs_tmp, "png16m", resolution=png_resolution)
            shutil.move(gs_tmp, png_path)
            png_ok = True
        except Exception as exc:
            if logger:
                logger.debug("Native PNG export failed for %s: %s", png_path,
                             exc)

    if not png_ok and pdf_ok and pdfium is None and logger:
        logger.debug(
            "pypdfium2 unavailable; PDF-to-PNG fallback unavailable.")

    # Python PDF fallback: if native PDF export failed but PNG exists,
    # build a PDF from the PNG with Pillow.
    if not pdf_ok and png_ok:
        if Image is None:
            try:
                from PIL import Image as pillow_image
                Image = pillow_image
            except Exception:
                Image = None

        if Image is not None:
            try:
                with Image.open(png_path) as im:
                    rgb = im.convert("RGB")
                    rgb.save(pdf_path, "PDF", resolution=float(png_resolution))
                pdf_ok = True
            except Exception as exc:
                if logger:
                    logger.debug("Pillow PNG-to-PDF fallback failed for %s: %s",
                                 pdf_path, exc)
        elif logger:
            logger.debug(
                "Pillow unavailable; PNG-to-PDF fallback unavailable.")

    return pdf_ok, png_ok


def export_text_rnamap_fallback(out_dir,
                                event_type,
                                map_name,
                                pdf_path,
                                png_path,
                                logger=None):
    """Create a simple Pillow plot from combined.RNAmap.txt when PyX export fails."""
    global Image

    if Image is None:
        try:
            from PIL import Image as pillow_image
            Image = pillow_image
        except Exception as exc:
            if logger:
                logger.debug("Pillow fallback unavailable: %s", exc)
            return False, False

    try:
        from PIL import ImageDraw, ImageFont
    except Exception as exc:
        if logger:
            logger.debug("Pillow drawing fallback unavailable: %s", exc)
        return False, False

    data_path = os.path.join(out_dir, "combined.RNAmap.txt")
    if not os.path.exists(data_path):
        if logger:
            logger.debug("RNA map text fallback input missing: %s", data_path)
        return False, False

    regions = []
    values = {}
    try:
        with open(data_path, newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                region = row.get("Region", "")
                if not region:
                    continue
                if region not in values:
                    values[region] = {"up": [], "down": [], "bg": []}
                    regions.append(region)
                values[region]["up"].append(float(row["AverageRawPeakCount_up"]))
                values[region]["down"].append(float(row["AverageRawPeakCount_down"]))
                values[region]["bg"].append(
                    float(row["AverageRawPeakCount_background"]))
    except Exception as exc:
        if logger:
            logger.debug("Could not read RNA map text fallback data: %s", exc)
        return False, False

    if not regions:
        return False, False

    width = 1500
    panel_h = 180
    top = 90
    bottom = 70
    left = 90
    right = 40
    gap = 34
    height = top + bottom + len(regions) * panel_h + (len(regions) - 1) * gap
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    font = ImageFont.load_default()

    colors = {
        "up": (205, 45, 45),
        "down": (45, 95, 205),
        "bg": (35, 35, 35),
    }

    title = f"{event_type} RNA map: {map_name}"
    draw.text((left, 24), title, fill=(20, 20, 20), font=font)
    legend_x = width - 355
    for idx, (name, color_value) in enumerate(
            [("Up", colors["up"]), ("Down", colors["down"]),
             ("Background", colors["bg"])]):
        y = 25 + idx * 20
        draw.line((legend_x, y + 7, legend_x + 32, y + 7),
                  fill=color_value,
                  width=3)
        draw.text((legend_x + 42, y), name, fill=(20, 20, 20), font=font)

    plot_w = width - left - right
    for region_index, region in enumerate(regions):
        y0 = top + region_index * (panel_h + gap)
        y1 = y0 + panel_h
        all_vals = (
            values[region]["up"] + values[region]["down"] +
            values[region]["bg"])
        max_val = max(all_vals) if all_vals else 1.0
        if max_val <= 0:
            max_val = 1.0

        draw.rectangle((left, y0, width - right, y1), outline=(210, 210, 210))
        draw.text((16, y0 + 6), region[:32], fill=(20, 20, 20), font=font)
        draw.text((left, y0 - 16), f"{max_val:.2f}", fill=(90, 90, 90),
                  font=font)
        draw.text((left, y1 + 2), "0", fill=(90, 90, 90), font=font)

        for series_name in ("bg", "down", "up"):
            series = values[region][series_name]
            if len(series) < 2:
                continue
            points = []
            denom = max(1, len(series) - 1)
            for idx, val in enumerate(series):
                x = left + int((idx / float(denom)) * plot_w)
                y = y1 - int((val / max_val) * (panel_h - 18)) - 9
                points.append((x, y))
            draw.line(points, fill=colors[series_name], width=2)

    png_ok = False
    pdf_ok = False
    try:
        image.save(png_path, format="PNG")
        png_ok = True
    except Exception as exc:
        if logger:
            logger.debug("Pillow RNA map PNG fallback failed: %s", exc)

    try:
        image.save(pdf_path, format="PDF")
        pdf_ok = True
    except Exception as exc:
        if logger:
            logger.debug("Pillow RNA map PDF fallback failed: %s", exc)

    return pdf_ok, png_ok


def export_motif_map_fallback(event_type,
                              map_name,
                              draw_points,
                              neg_pval_points,
                              max_point_value,
                              max_neg_pval,
                              pdf_path,
                              png_path,
                              logger=None):
    """Create a simple motif-map PDF/PNG when PyX font rendering is unavailable."""
    global Image

    if Image is None:
        try:
            from PIL import Image as pillow_image
            Image = pillow_image
        except Exception as exc:
            if logger:
                logger.debug("Pillow motif fallback unavailable: %s", exc)
            return False, False

    try:
        from PIL import ImageDraw, ImageFont
    except Exception as exc:
        if logger:
            logger.debug("Pillow motif drawing fallback unavailable: %s", exc)
        return False, False

    try:
        os.makedirs(os.path.dirname(png_path), exist_ok=True)
        os.makedirs(os.path.dirname(pdf_path), exist_ok=True)
    except Exception:
        pass

    region_count = len(draw_points)
    if region_count == 0:
        return False, False

    max_point_value = max(float(max_point_value or 0), 1e-12)
    max_neg_pval = max(float(max_neg_pval or 0), 1e-12)
    width = max(1100, region_count * 185 + 180)
    height = 520
    left = 70
    right = 35
    top = 95
    score_top = top
    score_bottom = 300
    pval_top = 345
    pval_bottom = 455
    gap = 10
    panel_width = max(60, int((width - left - right - gap * (region_count - 1)) / region_count))

    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    font = ImageFont.load_default()

    colors = {
        "up": (210, 45, 45),
        "down": (45, 85, 215),
        "bg": (40, 40, 40),
        "grid": (215, 215, 215),
    }

    safe_title = re.sub(r"\s+", " ", str(map_name))[:120]
    draw.text((left, 24), f"{event_type} motif map: {safe_title}",
              fill=(20, 20, 20), font=font)
    legend_items = [("Up", colors["up"]), ("Down", colors["down"]),
                    ("Background", colors["bg"]),
                    ("-log10 p up", colors["up"]),
                    ("-log10 p down", colors["down"])]
    legend_x = left
    for label, color_value in legend_items:
        y = 58
        draw.line((legend_x, y + 7, legend_x + 32, y + 7),
                  fill=color_value, width=3)
        draw.text((legend_x + 38, y), label, fill=(30, 30, 30), font=font)
        legend_x += 130

    draw.text((8, score_top + 60), "Motif score", fill=(80, 80, 80), font=font)
    draw.text((8, pval_top + 35), "-log10(p)", fill=(80, 80, 80), font=font)

    def y_from_value(value, max_value, y0, y1):
        frac = max(0.0, min(1.0, float(value) / max_value))
        return y1 - int(frac * (y1 - y0))

    def series_points(series, series_index, x0, x1, y0, y1, max_value):
        if len(series) < 2:
            return []
        denom = max(1, len(series) - 1)
        points = []
        for index, values in enumerate(series):
            x = x0 + int((index / float(denom)) * (x1 - x0))
            points.append((x, y_from_value(values[series_index], max_value, y0, y1)))
        return points

    for region_index, region_points in enumerate(draw_points):
        x0 = left + region_index * (panel_width + gap)
        x1 = x0 + panel_width
        for y0, y1 in ((score_top, score_bottom), (pval_top, pval_bottom)):
            draw.rectangle((x0, y0, x1, y1), outline=colors["grid"])
            draw.line((x0, (y0 + y1) // 2, x1, (y0 + y1) // 2),
                      fill=colors["grid"])
        draw.text((x0 + 4, score_bottom + 8), f"R{region_index + 1}",
                  fill=(80, 80, 80), font=font)

        for series_index, name in enumerate(("up", "down", "bg")):
            points = series_points(region_points, series_index, x0, x1,
                                   score_top + 8, score_bottom - 8,
                                   max_point_value)
            if len(points) >= 2:
                draw.line(points, fill=colors[name], width=2)

        if region_index < len(neg_pval_points):
            pval_region = neg_pval_points[region_index]
            for series_index, name in enumerate(("up", "down")):
                points = series_points(pval_region, series_index, x0, x1,
                                       pval_top + 8, pval_bottom - 8,
                                       max_neg_pval)
                if len(points) >= 2:
                    draw.line(points, fill=colors[name], width=1)

    png_ok = False
    pdf_ok = False
    try:
        image.save(png_path, format="PNG")
        png_ok = True
    except Exception as exc:
        if logger:
            logger.debug("Pillow motif PNG fallback failed: %s", exc)

    try:
        image.save(pdf_path, format="PDF")
        pdf_ok = True
    except Exception as exc:
        if logger:
            logger.debug("Pillow motif PDF fallback failed: %s", exc)

    return pdf_ok, png_ok


def suppress_pyx_text_cleanup(logger=None):
    """Avoid a second TeX timeout during PyX atexit cleanup after TeX failure."""
    try:
        engine = getattr(text, "defaulttextengine", None)
        instance = getattr(engine, "instance", None)
        cleanup = getattr(instance, "_cleanup", None)
        if cleanup is not None:
            atexit.unregister(cleanup)
            if logger:
                logger.debug("Suppressed PyX text cleanup after TeX failure")
            return True
    except Exception as exc:
        if logger:
            logger.debug("Could not suppress PyX text cleanup: %s", exc)
    return False
