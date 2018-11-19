#!/usr/bin/env python3
import sys 
import glob
import os

verbose = False

"""Add this script to your latexmkrc file:

example latexmkrc (assumes your images are in the figures folder and your main document is named 'figures':

$pdf_mode = 1;
$default_files = ('root.tex');
$bibtex_use = 2;
$clean_ext = 'thesis-blx.bib run.xml';

system('./topdf.py figures')
"""

def convert_to_pdf(folder):
    for svgname in glob.glob(os.path.join(folder,'*.svg')) + glob.glob(os.path.join(folder,'*/*.svg')):
        svgfile = os.path.basename(svgname)
        pdffile = ''.join(svgfile.split('.')[:-1])+'.pdf'
        pdfname = os.path.join(folder, pdffile)
        try:
            isModified = os.path.getmtime(svgname) > os.path.getmtime(pdfname)
        except FileNotFoundError:
            isModified = True
        if isModified:
            command = 'inkscape --export-pdf {} {}'.format(pdfname, svgname)
            print(command)
            os.system(command)
        elif verbose:
            print("{} is already converted and not modified, ignoring.".format(svgname))

if __name__ == '__main__':
    if len(sys.argv) < 2:
        raise Exception("You need to specify the image folder")
    if not os.path.exists(sys.argv[1]):
        raise Exception("Specified folder does not exist.")
    convert_to_pdf(sys.argv[1])
