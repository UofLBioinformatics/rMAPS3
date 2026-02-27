from pyx import *


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
