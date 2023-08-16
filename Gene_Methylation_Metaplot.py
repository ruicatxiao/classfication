#!/usr/bin/env python
import sys
import os
import time
from optparse import OptionParser
from collections import defaultdict

def main():
    parser = OptionParser()
    parser.add_option("-g", "--gff", dest="gff", help="gff input file, forced")
    parser.add_option("-a", "--allc", dest="allc", help="tsv input file, forced")
    parser.add_option("-o", dest="fOut", help="output directory, forced")
    (options, args) = parser.parse_args()
    
    if not (options.gff and options.allc and options.fOut):
        print("Error: Missing required arguments.")
        parser.print_help()
        sys.exit(1)
    
    if not os.path.exists(options.fOut):
        os.makedirs(options.fOut)
    
    winsize = 50
    binnum = 20
    context = {
        "CGA": "CG", "CGC": "CG", "CGG": "CG", "CGT": "CG", "CG": "CG", "CGN": "CG",
        "CHG": "CHG", "CHH": "CHH", "CAG": "CHG", "CCG": "CHG", "CTG": "CHG",
        "CAA": "CHH", "CAC": "CHH", "CAT": "CHH", "CCA": "CHH", "CCC": "CHH", "CCT": "CHH",
        "CTA": "CHH", "CTC": "CHH", "CTT": "CHH"
    }
    
    class_names = ["gene", "transposable_element_gene"]
    
    gene_bin_site_methylation(options.gff, options.allc, options.fOut, winsize, binnum, context, class_names)

def gene_bin_site_methylation(gff_file, allc_file, output_dir, winsize, binnum, context, class_names):
    up_max = winsize * binnum
    up_stop = binnum
    gb_start = binnum + 1
    gb_stop = binnum * 2
    gb_bin = binnum
    down_start = gb_stop + 1
    down_stop = down_start + binnum - 1
    chr_site = defaultdict(dict)
    chr = None
    
    with open(gff_file, 'r') as gff:
        for line in gff:
            line = line.strip()
            if line.startswith("#") or not line:
                continue
            columns = line.split("\t")
            if columns[0].startswith("Chr"):
                chr_match = re.match(r"Chr(\d+)", columns[0])
                if chr_match:
                    chr = chr_match.group(1)
                if chr and (columns[2] in class_names):
                    gene_name_match = re.search(r"Name=(\S+)", columns[-1])
                    if gene_name_match:
                        gene_name = gene_name_match.group(1)
                        for i in range(1, up_stop + 1):
                            binnum = i if columns[6] == "+" else down_stop + 1 - i
                            start = int(columns[3]) - up_max + (i - 1) * winsize
                            stop = start + winsize - 1
                            for j in range(start, stop + 1):
                                chr_site[(chr, j)] = (columns[2], binnum)
                        # More processing for gene body and downstream regions...
    
    mCbase = defaultdict(lambda: defaultdict(lambda: defaultdict(float)))
    allbase = defaultdict(lambda: defaultdict(lambda: defaultdict(float)))
    
    with open(allc_file, 'r') as allc:
        next(allc)  # Skip header
        for line in allc:
            line = line.strip()
            columns = line.split("\t")
            if (chr, int(columns[1])) in chr_site and context.get(columns[3]):
                group, binnum = chr_site[(chr, int(columns[1]))]
                mCbase[group][binnum][context[columns[3]]] += float(columns[4])
                allbase[group][binnum][context[columns[3]]] += float(columns[5])
    
    text = ["CG", "CHG", "CHH"]
    tsvname = os.path.splitext(os.path.basename(allc_file))[0]
    
    with open(os.path.join(output_dir, f"{tsvname}.bin.table"), 'w') as outfile:
        outfile.write("Group\tbinnum\tCG\tCHG\tCHH\n")
        for group in class_names:
            for binnum in sorted(allbase[group]):
                outfile.write(f"{group}\t{binnum}")
                for con in text:
                    if allbase[group][binnum][con]:
                        methylation_level = mCbase[group][binnum][con] / allbase[group][binnum][con]
                        outfile.write(f"\t{methylation_level}")
                    else:
                        outfile.write("\t0")
                outfile.write("\n")

if __name__ == "__main__":
    start_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    print(f"Program Starts Time: {start_time}")
    main()
    end_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    print(f"Program Ends Time: {end_time}\nDone. Total elapsed time: {time.time() - start_time}s")

