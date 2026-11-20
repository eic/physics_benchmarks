#!/usr/bin/env python3
"""Finalize one merged dihadron output: plots, histogram contents, benchmark tests.

    python3 finalize.py <merged.root> <plot dir> <hists.json> <tests.json> <label>

- <plot dir>/dihadron_<decade>.png: C(dphi) per trigger for each x_Bj decade that has
  triggers, all pairs and the CGC-weighted sum (gluon channels weighted with Cth).
- <hists.json>: for every histogram in the file (subdirectories included, keyed by path),
  its entries and bin contents over all cells, under/overflow included. A production
  system can digest these instead of keeping the ROOT file.
- <tests.json>: common_bench test format ("tests" list), read by bin/collect_tests.py.

The physics plots of the analysis (ep vs eA per decade, several samples) combine more
than one merged file and are not made here.
"""
import json
import math
import os
import sys

import ROOT

ROOT.gROOT.SetBatch(True)
# Set through the interpreter: assigning ROOT.gErrorIgnoreLevel from Python segfaults in
# the eic_xl 26.07.3 image (ROOT 6.38.04, Python 3.13).
ROOT.gInterpreter.ProcessLine('gErrorIgnoreLevel = kWarning;')
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

merged, plotdir, hists_out, tests_out, label = sys.argv[1:6]
os.makedirs(plotdir, exist_ok=True)

DECADES = {
    'x_1e-4_1e-3': '10^{-4} #leq #it{x} < 10^{-3}',
    'x_1e-3_1e-2': '10^{-3} #leq #it{x} < 10^{-2}',
    'x_1e-2_1e-1': '10^{-2} #leq #it{x} < 10^{-1}',
    'x_1e-1_1': '10^{-1} #leq #it{x} < 1',
}
Q = ['QCDC', 'LODIS', 'ResQQ', 'ResGQ']
G = ['PGF', 'ResQG', 'ResGG']
REBIN = 3

f = ROOT.TFile.Open(merged)
if not f or f.IsZombie():
    sys.exit(f'cannot open {merged}')


def walk(d, prefix=''):
    # Get(), not TKey.ReadObj(): ReadObj hands ownership of a subdirectory to Python, and
    # once Python frees it the file holds a dangling pointer (segfault on a later Get).
    for key in sorted(k.GetName() for k in d.GetListOfKeys()):
        o = d.Get(key)
        name = prefix + key
        if o.InheritsFrom('TDirectory'):
            yield from walk(o, name + '/')
        elif o.InheritsFrom('TH1'):
            yield name, o


hists = {}
for name, h in walk(f):
    hists[name] = {'entries': h.GetEntries(),
                   'bins': [h.GetBinContent(i) for i in range(h.GetNcells())]}
with open(hists_out, 'w') as out:
    json.dump(hists, out, separators=(',', ':'))

counters = f.Get('hCounters')
ntrig, events_read = counters.GetBinContent(1), counters.GetBinContent(4)

keep = []
plots = []
for dname, header in DECADES.items():
    d = f.Get(dname)
    nt = d.Get('hCounters').GetBinContent(1) if d else 0
    if not nt:
        continue
    allp = d.Get('deltaPhi').Clone(f'{dname}_all')
    sat = d.Get(f'deltaPhi{G[0]}Sat').Clone(f'{dname}_sat')
    for p in G[1:]:
        sat.Add(d.Get(f'deltaPhi{p}Sat'))
    for p in Q:
        sat.Add(d.Get(f'deltaPhi{p}'))
    c = ROOT.TCanvas(f'c_{dname}', '', 800, 600)
    c.SetLeftMargin(.14)
    ymax = 0
    for h, col, ms in ((allp, ROOT.kBlack, 20), (sat, ROOT.kRed + 1, 24)):
        h.SetDirectory(0)
        h.Rebin(REBIN)
        h.Scale(1. / nt / h.GetBinWidth(1))
        h.SetMarkerColor(col); h.SetLineColor(col); h.SetMarkerStyle(ms); h.SetMarkerSize(.8)
        ymax = max(ymax, h.GetMaximum())
        keep.append(h)
    allp.GetXaxis().SetTitle('#Delta#phi')
    allp.GetYaxis().SetTitle('C(#Delta#phi)'); allp.GetYaxis().SetTitleOffset(1.4)
    allp.SetMinimum(0); allp.SetMaximum(1.6 * ymax)
    allp.Draw('E1'); sat.Draw('E1 SAME')
    lg = ROOT.TLegend(.55, .75, .88, .88); lg.SetBorderSize(0); lg.SetFillStyle(0)
    lg.SetHeader(header)
    lg.AddEntry(allp, 'all pairs', 'lp'); lg.AddEntry(sat, 'Cth-weighted', 'lp')
    lg.Draw()
    # font 42 (size as a fraction of the pad): the repo's .rootlogon.C sets font 43, sized in pixels
    t = ROOT.TLatex(); t.SetNDC(); t.SetTextFont(42); t.SetTextSize(.033)
    keep += [lg, t, t.DrawLatex(.17, .85, label),
             t.DrawLatex(.17, .80, f'{int(nt)} triggers, {int(events_read)} events read')]
    png = os.path.join(plotdir, f'dihadron_{dname}.png')
    c.SaveAs(png)
    plots.append(png)


def test(name, title, description, quantity, target, value, ok):
    return {'name': name, 'title': title, 'description': description, 'quantity': quantity,
            'target': str(target), 'value': value, 'weight': 1.0,
            'result': 'pass' if ok else 'fail'}


tests = [
    test(f'{label}_events', 'Events read', 'The merged file counts the events it read',
         'events', '>0', events_read, events_read > 0),
    test(f'{label}_triggers', 'Triggers found', 'At least one trigger hadron in the merged file',
         'triggers', '>0', ntrig, ntrig > 0),
    test(f'{label}_plots', 'Plots made', 'At least one x decade has triggers and a plot',
         'plots', '>0', len(plots), len(plots) > 0),
]
with open(tests_out, 'w') as out:
    json.dump({'tests': tests}, out, indent=4)

print(f'finalize: {int(events_read)} events, {int(ntrig)} triggers, {len(hists)} histograms, '
      f'plots: {", ".join(os.path.basename(p) for p in plots)}')
sys.exit(0 if all(t['result'] == 'pass' for t in tests) else 1)
