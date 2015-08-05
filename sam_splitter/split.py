import argparse
import HTSeq
import os


def get_header(sam_file):
    header = []
    with open(sam_file, "r") as textfile:
        for line in textfile:
            if line.startswith("@"):
                header.append(line)
            else:
                break
    return header


def get_chroms(sam_file):
    sam_reader = HTSeq.SAM_Reader( sam_file )
    chroms = set()
    for a in sam_reader:
        if a.iv is not None:
            chroms.add(a.iv.chrom)
    return chroms


def new_file(chrom, header):
    a_file = open(chrom + ".sam", "w")
    for line in header:
        a_file.write(line)
    return a_file


def write_chroms(header, chroms, sam_file, report_file):
    files = {}
    counts = {}
    unaligned_file = new_file("unaligned", header)
    unaligned_count = 0
    for chrom in chroms:
        files[chrom] = new_file(chrom, header)
        counts[chrom] = 0
    sam_reader = HTSeq.SAM_Reader( sam_file )
    for a in sam_reader:
        if a.iv is not None:
            a_file = files[a.iv.chrom]
            a_file.write(a.get_sam_line())
            a_file.write("\n")
            counts[a.iv.chrom] += 1
        else:
            unaligned_file.write(a.get_sam_line())
            unaligned_file.write("\n")
            unaligned_count += 1
    with open(report_file, "w") as report:
        for chrom in sorted(chroms):
            report.write(str(chrom) + "  " + str(counts[chrom]) + "\n")
            files[chrom].close()
        report.write(str("unaligned") + "  " + str(unaligned_count) + "\n")
        unaligned_file.close()
    if unaligned_count == 0:
        os.remove("unaligned.sam")


def do_split(sam_file, report_file):
    header = get_header(sam_file)
    chroms = get_chroms(sam_file)
    write_chroms(header, chroms, sam_file, report_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Splits a sam file by chrom.')
    parser.add_argument("sam_file")
    parser.add_argument("report_file")

    # parse
    args = parser.parse_args()

    do_split(args.sam_file, args.report_file)
