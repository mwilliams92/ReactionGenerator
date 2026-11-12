#!/usr/bin/env python3
"""Extract the primary binary-emission residual excitation populations from TALYS.

The CSV contains a *probability mass per discrete level or continuum bin*, not
an energy-density. Populations are normalized to the sum of the extracted rows.
"""
import argparse
import csv
import math
from pathlib import Path
import re
import sys

SECTION_RE = re.compile(r'POPULATION\s+AFTER\s+BINARY\s+EMISSION', re.I)
NUCLIDE_RE = re.compile(
    r'Population\s+of\s+Z\s*=\s*(\d+)\s+N\s*=\s*(\d+)\s*'
    r'\(\s*([^)]*?)\s*\)\s+after\s+binary\s+(\S+)\s+emission\s*:\s*'
    r'([+\-\d.EeDd]+)', re.I)
META_RE = re.compile(r'Discrete\s+levels\s*:\s*(\d+)\s+Continuum\s+bins\s*:\s*(\d+)', re.I)
ROW_RE = re.compile(r'^\s*(\d+)\s+([+\-]?(?:\d+\.?\d*|\.\d+)(?:[EeDd][+\-]?\d+)?)\s+'
                    r'([+\-]?(?:\d+\.?\d*|\.\d+)(?:[EeDd][+\-]?\d+)?)(?=\s|$)')

def extract(text, z, a, ejectile, occurrence):
    lines = text.splitlines()
    sections = [i for i, line in enumerate(lines) if SECTION_RE.search(line)]
    # TALYS can be provided as a full output or as a standalone copied table.
    if not sections:
        sections = [-1]
    matches = []
    for si, start in enumerate(sections):
        stop = sections[si + 1] if si + 1 < len(sections) else len(lines)
        for i in range(start + 1, stop):
            m = NUCLIDE_RE.search(lines[i])
            if not m or (int(m[1]), int(m[1]) + int(m[2])) != (z, a):
                continue
            if m[4].lower() != ejectile.lower():
                continue
            meta = None
            rows = []
            begun = False
            for line in lines[i + 1:stop]:
                if NUCLIDE_RE.search(line) or line.lstrip().startswith(('++++++++', '########')):
                    break
                mm = META_RE.search(line)
                if mm:
                    meta = (int(mm[1]), int(mm[2]))
                r = ROW_RE.match(line)
                if r:
                    idx = int(r[1])
                    if rows and idx != rows[-1][0] + 1:
                        if begun:
                            break
                    rows.append((idx, float(r[2].replace('D','E').replace('d','e')),
                                 float(r[3].replace('D','E').replace('d','e'))))
                    begun = True
                elif begun and line.strip() and not line.lstrip().startswith(('J=', 'bin')):
                    break
                # Prevent accidentally ingesting the next section even if it starts with numbers.
                if meta and len(rows) >= meta[0] + 1 + meta[1]:
                    break
            if rows:
                matches.append((i + 1, m[3].strip(), float(m[5].replace('D','E').replace('d','e')), meta, rows))
    if not matches:
        raise ValueError(f'No binary {ejectile}-emission population found for Z={z}, A={a}.')
    if occurrence is None and len(matches) > 1:
        details = ', '.join(f'#{j}: line {m[0]}' for j, m in enumerate(matches, 1))
        raise ValueError(f'Found {len(matches)} matching tables ({details}). '
                         'Specify --occurrence N for the intended TALYS incident energy.')
    which = (occurrence or 1) - 1
    if which < 0 or which >= len(matches):
        raise ValueError(f'Occurrence {occurrence} not available ({len(matches)} matching tables).')
    return matches[which]

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument('talys_output', type=Path, help='TALYS standard output text file')
    ap.add_argument('--Z', '-Z', type=int, required=True, help='Residual proton number')
    ap.add_argument('--A', '-A', type=int, required=True, help='Residual mass number')
    ap.add_argument('--ejectile', default='neutron', help='Binary ejectile name in TALYS table')
    ap.add_argument('-o', '--output', type=Path, help='Output CSV path')
    ap.add_argument('--occurrence', type=int, help='Select matching table (1-based), for multi-energy output')
    ap.add_argument('--discrete-last-bin', type=int, help='Override last discrete bin index from TALYS header')
    args = ap.parse_args()
    if args.A < args.Z or args.Z < 0:
        ap.error('Require 0 <= Z <= A')
    try:
        line_no, isotope, header_total, meta, rows = extract(
            args.talys_output.read_text(errors='replace'), args.Z, args.A,
            args.ejectile, args.occurrence)
        if not meta and args.discrete_last_bin is None:
            raise ValueError('Could not find "Discrete levels: ... Continuum bins: ..." header; '
                             'use --discrete-last-bin to specify boundary.')
        last_discrete = (args.discrete_last_bin if args.discrete_last_bin is not None else meta[0])
        total = sum(p for _, _, p in rows)
        if total <= 0 or not math.isfinite(total) or any(p < 0 or not math.isfinite(p) for _,_,p in rows):
            raise ValueError('Populations are nonfinite, negative, or sum to zero.')
        outfile = args.output or Path(f'population_Z{args.Z}_A{args.A}.csv')
        with outfile.open('w', newline='') as fh:
            writer = csv.writer(fh)
            writer.writerow(['bin', 'type', 'Ex_MeV', 'population_mb', 'probability', 'cumulative_probability'])
            cumulative = 0.0
            for bin_id, ex, pop in rows:
                cumulative += pop / total
                writer.writerow([bin_id, 'discrete' if bin_id <= last_discrete else 'continuum',
                                 f'{ex:.8g}', f'{pop:.10g}', f'{pop/total:.12g}', f'{min(cumulative, 1.0):.12g}'])
        print(f'{isotope} (Z={args.Z}, A={args.A}), binary {args.ejectile} emission')
        print(f'Matching table at line {line_no}; {len(rows)} rows; last discrete bin {last_discrete}')
        if meta and len(rows) != meta[0] + 1 + meta[1]:
            print(f'WARNING: Expected {meta[0] + 1 + meta[1]} rows, found {len(rows)}', file=sys.stderr)
        print(f'Summed population: {total:.8g} mb; TALYS heading: {header_total:.8g} mb')
        if header_total > 0 and abs(total-header_total)/header_total > 0.01:
            print('WARNING: Row sum differs from table heading by >1%; inspect the TALYS table.', file=sys.stderr)
        print(f'CSV: {outfile}')
    except (OSError, ValueError) as exc:
        ap.exit(1, f'Error: {exc}\n')

if __name__ == '__main__':
    main()
