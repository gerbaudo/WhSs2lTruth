#!/bin/env python

# plot, in the (mc1,mn1) plane the signal acceptance and efficiency
#
# acceptance = N_{fiducial} / N_{generated}
#            = (N events in SR with selection on truth objects) / ((N_{events generated} / (BR * filter_eff))
# efficiency = N_{fiducial-reco} / N_{fiducial}
#            = (N events in SR with selection on reco objects) / (N events in SR with selection on truth objects)
#
# For more details on these definitions, see
# https://indico.cern.ch/event/159188/material/slides/0?contribId=37
# and
# http://arxiv.org/pdf/1403.4853v1.pdf (appendix A)
#
# Nov  6, first version
# Nov 10, update with numbers from Josie
# Nov 12, update: using TGraph2D
# Dec 09, update: grey palette, adjust font sizes, labels
# Dec 12, bugfix: missing parentheses acceptance denominator
#
#___________________________________________________________

import array
import ROOT as r
r.gROOT.SetBatch(True)                     # no windows popping up
r.PyConfig.IgnoreCommandLineOptions = True # don't let root steal our cmd-line options

def setAtlasStyle() :
    aStyle = getAtlasStyle()
    r.gROOT.SetStyle("ATLAS")
    r.gROOT.ForceStyle()
    r.gStyle.SetPalette(52)
    nElementsPalette = 40
    r.gStyle.SetPalette(nElementsPalette, array.array('i', [61+i for i in range(nElementsPalette)]))

def getAtlasStyle() :
    style = r.TStyle('ATLAS', 'Atlas style')
    white = 0
    style.SetFrameBorderMode(white)
    style.SetFrameFillColor(white)
    style.SetCanvasBorderMode(white)
    style.SetCanvasColor(white)
    style.SetPadBorderMode(white)
    style.SetPadColor(white)
    style.SetStatColor(white)
    #style.SetPaperSize(20,26)
    style.SetPadTopMargin(0.05)
    style.SetPadRightMargin(0.05)
    nokidding = 0.75 # limit the exaggerated margins
    style.SetPadBottomMargin(nokidding*0.16)
    style.SetPadLeftMargin(nokidding*0.16)
    style.SetPadRightMargin(nokidding*0.10)
    style.SetTitleXOffset(nokidding*1.25)
    style.SetTitleYOffset(nokidding*1.25)
    style.SetTitleXSize(1.25*style.GetTitleSize())
    style.SetTitleYSize(1.25*style.GetTitleSize())
    font, fontSize = 42, 0.04 # helvetica large
    style.SetTextFont(font)
    #style.SetTextSize(fontSize)
    style.SetLabelFont(font,"xyz")
    style.SetTitleFont(font,"xyz")
    style.SetPadTickX(1)
    style.SetPadTickY(1)
    style.SetOptStat(0)
    style.SetOptTitle(0)
    style.SetEndErrorSize(0)
    return style

def topRightLabel(pad, label, xpos=None, ypos=None, align=33) :
    pad.cd()
    tex = r.TLatex(0.0, 0.0, '')
    tex.SetTextFont(42)
    tex.SetTextSize(0.75*tex.GetTextSize())
    tex.SetNDC()
    tex.SetTextAlign(align)
    tex.DrawLatex((1.0-1.1*pad.GetRightMargin()) if not xpos else xpos,
                  (1.0-1.1*pad.GetTopMargin()) if not ypos else ypos,
                  label)
    if hasattr(pad, '_labels') : pad._labels.append(tex)
    else : pad._labels = [tex]
    return tex

def topLeftLabel(pad, label, xpos=None, ypos=None, align=13) :
    pad.cd()
    tex = r.TLatex(0.0, 0.0, '')
    tex.SetTextFont(42)
    tex.SetTextSize(0.75*tex.GetTextSize())
    tex.SetNDC()
    tex.SetTextAlign(align)
    tex.DrawLatex((0.0+1.1*pad.GetLeftMargin()) if not xpos else xpos,
                  (1.0-1.1*pad.GetTopMargin()) if not ypos else ypos,
                  label)
    if hasattr(pad, '_labels') : pad._labels.append(tex)
    else : pad._labels = [tex]
    return tex

def drawAtlasLabel(pad, xpos=None, ypos=None) :
    label = "#bf{#it{ATLAS}} Internal"
    return topLeftLabel(pad, label, xpos, ypos)

def nicelabel(l):
    return {'sr1jmm':'SR#mu#mu-1',
            'sr1jem':'SRe#mu-1',
            'sr1jee':'SRee-1',
            'sr2jmm':'SR#mu#mu-2',
            'sr2jem':'SRe#mu-2',
            'sr2jee':'SRee-2',
        }[l]

# main
#___________________________________________________________

counts = [
    # from https://svnweb.cern.ch/trac/atlasoff/browser/PhysicsAnalysis/SUSYPhys/SUSYTools/trunk/data/mc12_8TeV/Herwigpp_UEEE3_CTEQ6L1_simplifiedModel_wA.txt
    # N_generated from ami
    # N_fiducial, N_fiducial-reco from Josie, see email 2014-11-10
    # dsID    MC1,MN2[GeV]   MN1[GeV]   XS    BF                 filterEff      N_generated    {N_fiducial  N_fiducial-reco} x sr
    (177501,  130.0,    0.0,  4.2439701557,   0.3063600000,  0.1587900000,     99999,  {'sr1jee':(429.0, 122.0), 'sr1jmm':(433.0,  111.0), 'sr1jem':(608.0,  319.0), 'sr2jee':(207.0,  105.0), 'sr2jmm':(160.0,  98.0 ), 'sr2jem':(283.0,  229.0)}),
    (177502,  140.0,   10.0,  3.1979225874,   0.3063600000,  0.1606000000,     80000,  {'sr1jee':(369.0, 118.0), 'sr1jmm':(379.0,  102.0), 'sr1jem':(472.0,  257.0), 'sr2jee':(181.0,  83.0 ), 'sr2jmm':(138.0,  71.0 ), 'sr2jem':(248.0,  187.0)}),
    (177503,  150.0,    0.0,  2.4526928067,   0.3063600000,  0.1625300000,    100000,  {'sr1jee':(490.0, 131.0), 'sr1jmm':(514.0,  140.0), 'sr1jem':(609.0,  358.0), 'sr2jee':(270.0,  108.0), 'sr2jmm':(171.0,  110.0), 'sr2jem':(355.0,  244.0)}),
    (177504,  152.5,   22.5,  2.3024949431,   0.3063600000,  0.1632400000,     80000,  {'sr1jee':(369.0, 106.0), 'sr1jmm':(410.0,  107.0), 'sr1jem':(540.0,  246.0), 'sr2jee':(205.0,  85.0 ), 'sr2jmm':(140.0,  84.0 ), 'sr2jem':(277.0,  199.0)}),
    (177505,  162.5,   12.5,  1.8069593906,   0.3063600000,  0.1647800000,     79999,  {'sr1jee':(399.0, 123.0), 'sr1jmm':(381.0,  110.0), 'sr1jem':(508.0,  273.0), 'sr2jee':(221.0,  103.0), 'sr2jmm':(132.0,  77.0 ), 'sr2jem':(271.0,  226.0)}),
    (177506,  175.0,    0.0,  1.3526386917,   0.3063600000,  0.1680500000,    100000,  {'sr1jee':(545.0, 157.0), 'sr1jmm':(467.0,  116.0), 'sr1jem':(642.0,  426.0), 'sr2jee':(284.0,  143.0), 'sr2jmm':(164.0,  92.0 ), 'sr2jem':(350.0,  295.0)}),
    (177507,  165.0,   35.0,  1.7067664266,   0.3063600000,  0.1639400000,     80000,  {'sr1jee':(415.0, 87.0 ), 'sr1jmm':(410.0,  102.0), 'sr1jem':(523.0,  249.0), 'sr2jee':(206.0,  98.0 ), 'sr2jmm':(157.0,  86.0 ), 'sr2jem':(289.0,  196.0)}),
    (177508,  175.0,   25.0,  1.3526386917,   0.3063600000,  0.1663000000,     80000,  {'sr1jee':(415.0, 130.0), 'sr1jmm':(390.0,  112.0), 'sr1jem':(517.0,  314.0), 'sr2jee':(234.0,  120.0), 'sr2jmm':(128.0,  81.0 ), 'sr2jem':(287.0,  228.0)}),
    (177509,  187.5,   12.5,  1.0351453722,   0.3063600000,  0.1682700000,     79999,  {'sr1jee':(447.0, 159.0), 'sr1jmm':(399.0,  124.0), 'sr1jem':(506.0,  359.0), 'sr2jee':(242.0,  128.0), 'sr2jmm':(141.0,  86.0 ), 'sr2jem':(312.0,  243.0)}),
    (177510,  200.0,    0.0,  0.8022477627,   0.3063600000,  0.1705300000,    100000,  {'sr1jee':(523.0, 209.0), 'sr1jmm':(496.0,  178.0), 'sr1jem':(626.0,  477.0), 'sr2jee':(351.0,  176.0), 'sr2jmm':(181.0,  108.0), 'sr2jem':(359.0,  327.0)}),
    (177511,  177.5,   47.5,  1.2798982263,   0.3063600000,  0.1666300000,     49999,  {'sr1jee':(250.0, 72.0 ), 'sr1jmm':(246.0,  64.0 ), 'sr1jem':(331.0,  184.0), 'sr2jee':(132.0,  59.0 ), 'sr2jmm':(88.0 ,  68.0 ), 'sr2jem':(154.0,  138.0)}),
    (177512,  187.5,   37.5,  1.0351453722,   0.3063600000,  0.1674600000,     49999,  {'sr1jee':(268.0, 93.0 ), 'sr1jmm':(261.0,  85.0 ), 'sr1jem':(302.0,  235.0), 'sr2jee':(129.0,  85.0 ), 'sr2jmm':(79.0 ,  49.0 ), 'sr2jem':(195.0,  144.0)}),
    (177513,  200.0,   25.0,  0.8022477627,   0.3063600000,  0.1704700000,     50000,  {'sr1jee':(253.0, 90.0 ), 'sr1jmm':(226.0,  71.0 ), 'sr1jem':(326.0,  252.0), 'sr2jee':(173.0,  61.0 ), 'sr2jmm':(78.0 ,  50.0 ), 'sr2jem':(170.0,  136.0)}),
    (177514,  212.5,   12.5,  0.6294424236,   0.3063600000,  0.1729700000,     50000,  {'sr1jee':(267.0, 120.0), 'sr1jmm':(250.0,  81.0 ), 'sr1jem':(341.0,  266.0), 'sr2jee':(173.0,  100.0), 'sr2jmm':(97.0 ,  52.0 ), 'sr2jem':(215.0,  181.0)}),
    (177515,  225.0,    0.0,  0.5012571812,   0.3063600000,  0.1746500000,    100000,  {'sr1jee':(569.0, 233.0), 'sr1jmm':(476.0,  177.0), 'sr1jem':(655.0,  515.0), 'sr2jee':(354.0,  183.0), 'sr2jmm':(188.0,  114.0), 'sr2jem':(380.0,  347.0)}),
    (177516,  190.0,   60.0,  0.9834355116,   0.3063600000,  0.1684600000,     50000,  {'sr1jee':(235.0, 73.0 ), 'sr1jmm':(265.0,  83.0 ), 'sr1jem':(349.0,  167.0), 'sr2jee':(144.0,  67.0 ), 'sr2jmm':(101.0,  46.0 ), 'sr2jem':(169.0,  127.0)}),
    (177517,  200.0,   50.0,  0.8022477627,   0.3063600000,  0.1683700000,     50000,  {'sr1jee':(253.0, 93.0 ), 'sr1jmm':(241.0,  77.0 ), 'sr1jem':(336.0,  216.0), 'sr2jee':(152.0,  74.0 ), 'sr2jmm':(82.0 ,  57.0 ), 'sr2jem':(184.0,  155.0)}),
    (177518,  212.5,   37.5,  0.6294424236,   0.3063600000,  0.1730100000,     49900,  {'sr1jee':(293.0, 94.0 ), 'sr1jmm':(257.0,  94.0 ), 'sr1jem':(314.0,  219.0), 'sr2jee':(155.0,  72.0 ), 'sr2jmm':(106.0,  64.0 ), 'sr2jem':(188.0,  183.0)}),
    (177519,  225.0,   25.0,  0.5012571812,   0.3063600000,  0.1752300000,     50000,  {'sr1jee':(263.0, 112.0), 'sr1jmm':(255.0,  96.0 ), 'sr1jem':(370.0,  264.0), 'sr2jee':(185.0,  100.0), 'sr2jmm':(100.0,  56.0 ), 'sr2jem':(216.0,  191.0)}),
    (177520,  237.5,   12.5,  0.4004180953,   0.3063600000,  0.1764100000,     50000,  {'sr1jee':(271.0, 132.0), 'sr1jmm':(217.0,  99.0 ), 'sr1jem':(336.0,  267.0), 'sr2jee':(173.0,  101.0), 'sr2jmm':(108.0,  69.0 ), 'sr2jem':(207.0,  195.0)}),
    (177521,  250.0,    0.0,  0.3232687637,   0.3063600000,  0.1783900000,     99900,  {'sr1jee':(622.0, 244.0), 'sr1jmm':(479.0,  187.0), 'sr1jem':(676.0,  511.0), 'sr2jee':(370.0,  239.0), 'sr2jmm':(195.0,  133.0), 'sr2jem':(386.0,  386.0)}), # was 390, buggy?
    (177522,  202.5,   72.5,  0.7630156279,   0.3063600000,  0.1680100000,     49900,  {'sr1jee':(247.0, 95.0 ), 'sr1jmm':(254.0,  70.0 ), 'sr1jem':(328.0,  198.0), 'sr2jee':(153.0,  73.0 ), 'sr2jmm':(79.0 ,  51.0 ), 'sr2jem':(148.0,  131.0)}),
    (177523,  212.5,   62.5,  0.6294424236,   0.3063600000,  0.1698000000,     50000,  {'sr1jee':(289.0, 88.0 ), 'sr1jmm':(261.0,  77.0 ), 'sr1jem':(310.0,  202.0), 'sr2jee':(152.0,  82.0 ), 'sr2jmm':(103.0,  66.0 ), 'sr2jem':(188.0,  160.0)}),
    (177524,  225.0,   50.0,  0.5012571812,   0.3063600000,  0.1738200000,     50000,  {'sr1jee':(280.0, 111.0), 'sr1jmm':(235.0,  74.0 ), 'sr1jem':(346.0,  226.0), 'sr2jee':(169.0,  91.0 ), 'sr2jmm':(90.0 ,  52.0 ), 'sr2jem':(221.0,  156.0)}),
    (177525,  237.5,   37.5,  0.4004180953,   0.3063600000,  0.1773600000,     50000,  {'sr1jee':(282.0, 132.0), 'sr1jmm':(254.0,  95.0 ), 'sr1jem':(336.0,  258.0), 'sr2jee':(161.0,  97.0 ), 'sr2jmm':(94.0 ,  54.0 ), 'sr2jem':(208.0,  201.0)}),
    (177526,  250.0,   25.0,  0.3232687637,   0.3063600000,  0.1785500000,     50000,  {'sr1jee':(332.0, 118.0), 'sr1jmm':(247.0,  103.0), 'sr1jem':(347.0,  278.0), 'sr2jee':(177.0,  106.0), 'sr2jmm':(90.0 ,  43.0 ), 'sr2jem':(220.0,  194.0)}),
    (177527,  262.5,   12.5,  0.2637965977,   0.3063600000,  0.1789700000,     50000,  {'sr1jee':(309.0, 128.0), 'sr1jmm':(236.0,  94.0 ), 'sr1jem':(324.0,  263.0), 'sr2jee':(188.0,  106.0), 'sr2jmm':(89.0 ,  49.0 ), 'sr2jem':(215.0,  202.0)}),
# missing    (177528,      275.0,         0.0,    0.2168904692,   0.3063600000,    0.1821800000,         95000,           1,        1),
]

def dsid_cal(line) : return line[0]
def mc1_val(line) : return line[1]
def mn1_val(line) : return line[2]
def xsec_val(line) : return line[3]
def bf_val(line): return line[4]
def filter_eff_val(line) : return line[5]
def n_gen_val(line) : return line[6]
def n_fiducial_val(line, region) : return line[7][region][0]
def n_fiducial_reco_val(line, region) : return line[7][region][1]


mc1Range = {'min': min([mc1_val(c) for c in counts]), 'max' : max([mc1_val(c) for c in counts])}
mn1Range = {'min': min([mn1_val(c) for c in counts]), 'max' : max([mn1_val(c) for c in counts])}
mc1Range = {'max': 270.0, 'min': 130.0} # set the range by hand (auto is ok as a first approx)
mn1Range = {'max': 80.0, 'min': 0.0}

setAtlasStyle()

for selection in ['sr1jee', 'sr1jmm', 'sr1jem', 'sr2jee', 'sr2jmm', 'sr2jem']:
    title  = ''
    title += '; m_{#tilde{#chi}_{1}^{#pm},#tilde{#chi}_{2}^{0}} [GeV]'
    title += '; m_{#tilde{#chi}_{1}^{0}} [GeV]'
    
    histo_pad_master = r.TH2F('acceptance_'+selection, title,
                              100, float(mc1Range['min']), float(mc1Range['max']),
                              100, float(mn1Range['min']), float(mn1Range['max']))
    histo_acceptance = r.TH2F('acceptance_'+selection, title,
                              100, float(mc1Range['min']), float(mc1Range['max']),
                              100, float(mn1Range['min']), float(mn1Range['max']))
    histo_efficiency = r.TH2F('efficiency_'+selection, title,
                              100, float(mc1Range['min']), float(mc1Range['max']),
                              100, float(mn1Range['min']), float(mn1Range['max']))
    percent = 100.
    acceptance_scale_factor = 1.0e4
    acceptance_scale_label = '10^{4}'
    for count in counts:
        mc1, mn1 = mc1_val(count), mn1_val(count)
        bf, filterEff = bf_val(count), filter_eff_val(count)
        N_generated = n_gen_val(count)
        N_fiducial, N_fiducial_reco = n_fiducial_val(count, selection), n_fiducial_reco_val(count, selection)

        histo_acceptance.Fill(mc1, mn1, N_fiducial / (N_generated / (bf * filterEff)) * acceptance_scale_factor)
        histo_efficiency.Fill(mc1, mn1, N_fiducial_reco / N_fiducial * percent)
    # plot
    r.gStyle.SetPaintTextFormat('.3f')
    maxEff = 100.
    for h in [histo_acceptance, histo_efficiency]:
        c = r.TCanvas('c_'+h.GetName(), '', 800, 600)
        c.cd()
        c.SetRightMargin(2.0*c.GetRightMargin())
        graph = r.TGraph2D(h)
        h.SetStats(0)
        h.SetMarkerSize(1.5*h.GetMarkerSize())
        # histo_pad_master.SetMaximum(h.GetMaximum())
        histo_pad_master.Draw('axis') # just to get axes and palette
	graph.Draw("colz same")
 	c.Update()

	text_eff = []
	text_acc = []
	counter = 0
	for count in counts:
	  mc1_txt, mn1_txt = mc1_val(count), mn1_val(count)
	  bf_txt, filterEff_txt = bf_val(count), filter_eff_val(count)
	  N_generated_txt = n_gen_val(count)
	  N_fiducial_txt, N_fiducial_reco_txt = n_fiducial_val(count, selection), n_fiducial_reco_val(count, selection)
	  acceptance_txt = N_fiducial_txt / (N_generated_txt / (bf_txt * filterEff_txt)) * acceptance_scale_factor
	  efficiency_txt = N_fiducial_reco_txt / N_fiducial_txt * percent
	  #print N_fiducial_reco_txt, N_fiducial_txt
	  if "accep" in c.GetName():
	    text_acc.append(r.TText(mc1_txt, mn1_txt, str(round(acceptance_txt,1))))
	    text_acc[counter].SetTextSize(0.03)
	    # text_acc[counter].SetTextAngle(15)
	    text_acc[counter].Draw("same")
	    #print mc1_txt, mn1_txt, acceptance_txt, str(round(acceptance_txt,2))
	  if "efficiency" in c.GetName():
	    text_eff.append(r.TText(mc1_txt, mn1_txt, str(int(round(efficiency_txt,0)))))
	    text_eff[counter].SetTextSize(0.03)
	    # text_eff[counter].SetTextAngle(15)
	    text_eff[counter].Draw("same")
	    #print mc1_txt, mn1_txt, efficiency_txt, str(round(efficiency_txt,2))
	  counter += 1

	c.Update()
	xAx = histo_pad_master.GetXaxis()
        yAx = histo_pad_master.GetYaxis()
	xAx.SetTitle("m_{#tilde{#chi}_{1}^{#pm},#tilde{#chi}_{2}^{0}} [GeV]")
	yAx.SetTitle("m_{#tilde{#chi}_{1}^{0}} [GeV]")
	# xAx.SetLimits(125, 350)
	# yAx.SetLimits(0, 90)
	xAx.SetRangeUser(125.0, 350.0)
	yAx.SetRangeUser(0.0, 90.0)
        minTitleSize = min(a.GetTitleSize() for a in [xAx, yAx])
        xAx.SetTitleSize(minTitleSize)
        yAx.SetTitleSize(minTitleSize)
	xAx.SetTitleOffset(1.1*xAx.GetTitleOffset())
	yAx.SetTitleOffset(1.1*yAx.GetTitleOffset())
        # zAx.SetTitle("Acceptance [%]" if 'accep' in c.GetName() else 'Efficiency [%]')
	c.Update()

        zAx = graph.GetZaxis()
        zAx.SetTitle((acceptance_scale_label+' x Acceptance') if 'accep' in c.GetName() else 'Efficiency [%]')
        zAx.SetTitleOffset(1.1*zAx.GetTitleOffset())

        # labels etc.
        drawAtlasLabel(c, ypos=0.92)
        topLeftLabel(c, '#sqrt{s} = 8 TeV, 20.3 fb^{-1}', ypos=0.86)
        topRightLabel(c, nicelabel(selection), xpos=0.75, ypos=0.92)
        c.Update()

        c.SaveAs('../pics_acc_eff/' + c.GetName()+'.eps')
        c.SaveAs('../pics_acc_eff/' + c.GetName()+'.png')
