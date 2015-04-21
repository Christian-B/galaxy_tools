from collections import Counter, OrderedDict
from optparse import OptionParser, OptionGroup
import re
import sys
     
#Based on https://samtools.github.io/hts-specs/SAMv1.pdf     
     
def int_value(line_number, line, tag, value, min_value, max_value):
    try:
        result = int(value)
    except:
        print >> sys.stderr, "Line", line_number, "has an none ingeter", tag, "value"
        print >> sys.stderr, "Found",line
        sys.exit(1)
    if result < min_value:
        print >> sys.stderr, "Line", line_number, "has a", tag, "value below", min_value
        print >> sys.stderr, "Found",line
        sys.exit(1)
    if result > max_value:
        print >> sys.stderr, "Line", line_number, "has a", tag, "value above", max_value
        print >> sys.stderr, "Found",line
        sys.exit(1)
    return result        
   
sequence_dictionary = OrderedDict()     
   
class Sequence:
    def __init__(self, line_number, line):
        global _counterKeys
        self.sequence_counter = Counter()
        self.quality_counter = Counter()
        self.cigar_counter = Counter()
        self.length = None
        self.count = 0
        self.max_mapq = None
        self.min_mapq = None
        if line_number:    
            self.name = None
            parts = line.strip().split("\t")
            for part in parts[1:]:
                tag = part[:3]
                value = part[3:]
                if tag == "SN:":
                    self.name = value
                elif tag == "LN:":
                    self.length = int_value(line_number, line, "LN:", value, 1, 2**31-1)           
                elif tag == "AS:":
                    self.assembly = value             
                elif tag == "M5:":
                    self.checksum = value             
                elif tag == "SP:":
                    self.species = value             
                elif tag == "UR:":
                    self.uri = value             
                else:
                    print >> sys.stderr, "Line", line_number, "has an unexpected tag",tag
                    print >> sys.stderr, "Found",line
                    sys.exit(1)
            if not self.name:    
                print >> sys.stderr, "Line", line_number, "has a sequence('@SQ') with no name('SN:')"
                print >> sys.stderr, "Found",line
                sys.exit(1)
            if not self.length:    
                print >> sys.stderr, "Line", line_number, "has a sequence('@SQ') with no length('LN:')"
                print >> sys.stderr, "Found",line
                sys.exit(1)    
        else:
            self.name = line
            _sequence_counter_keys = None
            
    def addAlignment(self, alignment):
        self.count += 1
        self.sequence_counter.update(alignment.seq)
        self.quality_counter.update(alignment.qual)
        self.cigar_counter.update(alignment.cigar_long)
        if self.max_mapq:
            if self.max_mapq < alignment.mapq:
                self.max_mapq = alignment.mapq       
            if self.min_mapq > alignment.mapq:
                self.min_mapq = alignment.mapq       
        else:
            self.max_mapq = alignment.mapq        
            self.min_mapq = alignment.mapq       
        if self != _all_sequence:
            _all_sequence.addAlignment(alignment)               
        else:
            _sequence_keys = None

    def __str__(self):
        reply = str(self.count) 
        reply+=  "\t" 
        reply+= str(self.min_mapq) 
        reply+= "\t" 
        reply+= str(self.max_mapq)
        for key in sequenceKeys():            
            reply+= "\t"
            reply+= str(self.sequence_counter[key]) 
        reply+= "\t" 
        reply+= str(sum(self.sequence_counter.values()) )  
        for key in qualityKeys():            
            reply+= "\t"
            reply+= str(self.quality_counter[key]) 
        reply+= "\t" 
        reply+= str(sum(self.quality_counter.values()))   
        for key in cigarKeys():            
            reply+= "\t"
            reply+= str(self.cigar_counter[key]) 
        reply+= "\t" 
        reply+= str(sum(self.cigar_counter.values()))   
        reply+= "\t" 
        reply+= self.name 
        return reply      

_all_sequence =  Sequence(None, "All")  
sequence_dictionary["*"] = Sequence(None, "Unkown")           

_sequence_keys  = None                         
_quality_keys  = None                         
_cigar_keys  = None                         

def updateKeys():
    global _sequence_keys, _quality_keys, _cigar_keys
    _sequence_keys = sorted(_all_sequence.sequence_counter.keys())
    _quality_keys = sorted(_all_sequence.quality_counter.keys())
    _cigar_keys = sorted(_all_sequence.cigar_counter.keys())
    
def sequenceKeys():
    if not _sequence_keys:
        updateKeys()
    return _sequence_keys
            
def qualityKeys():
    if not _sequence_keys:
        updateKeys()
    return _quality_keys
            
def cigarKeys():
    if not _sequence_keys:
        updateKeys()
    return _cigar_keys
            
def strHeader():
    line = ["Count","Min Q","max Q"]
    line.extend(sequenceKeys())
    line.append("TSeq")
    line.extend(qualityKeys())
    line.append("TQual")
    line.extend(cigarKeys())
    line.append("Tcigar")
    line.append("Name")
    return "\t".join(line)  

program_dictionary = OrderedDict()     

class Program:
    def __init__(self, line_number, line):
        self.id = None
        parts = line.strip().split("\t")
        for part in parts[1:]:
            tag = part[:3]
            value = part[3:]
            if tag == "ID:":
                self.id = value
            elif tag == "PN:":
                self.name = value
            elif tag == "CL:":
                self.command_line = value             
            elif tag == "PP:":
                if program_dictionary[value]:
                    self.previous = program_dictionary[value] 
                else:
                    print >> sys.stderr, "Program(@PG) line", line_number, "has an unused PP id:",value
                    print >> sys.stderr, "Found",line
                    sys.exit(1)
            elif tag == "DS:":
                self.description = value             
            elif tag == "VN:":
                self.version = value             
            else:
                print >> sys.stderr, "Program(@PG) line", line_number, "has an unexpected tag",tag
                print >> sys.stderr, "Found",line
                sys.exit(1)
        if not self.id:    
            print >> sys.stderr, "Line(@PG)", line_number, "has no id('PG:')"
            print >> sys.stderr, "Found",line
            sys.exit(1)
                 
class Alignment():
    def __init__(self, line_number, line):
        parts = line.strip().split("\t")
        if len(parts) < 11:
            print >> sys.stderr, "Line", line_number, "has an allignnment with less than 11 columns."
            print >> sys.stderr, "Found",line
            sys.exit(1)
        self.qname = parts[0]    
      
        flag = int_value(line_number, line, "FLAG", parts[1], 0, 2**16-1)
        self.multiple_segments = flag & 1 != 0
        self.properly_aligned = flag & 2 != 0
        self.unmapped = flag & 4 != 0
        self.next_unmapped = flag & 4 != 0
        self.reverse_complemented = flag & 16 != 0  #0x10
        self.next_reverse_complemented = flag & 32 != 0 #OX20
        self.first_segment = flag & 64 != 0 #Ox40
        self.last_segment = flag & 128 != 0 #Ox80
        self.secondary = flag & 256 != 0 #Ox100
        self.failed_quality_control = flag & 512 != 0 #Ox200
        self.duplicate = flag & 1024 != 0 #Ox400
        self.supplement = flag & 2048 != 0 #Ox800
      
        self.rname = parts[2]              
        self.ps = int_value(line_number, line, "POS", parts[3], 0, 2**31-1)
        self.mapq = int_value(line_number, line, "MAPQ", parts[4], 0, 2**8-1)
       
        self.cigar = parts[5]
        if self.cigar != "*":
            cigars = re.findall("[0-9]+[MIDNSHPX=]",self.cigar)
            if len(self.cigar) != sum(len(c) for c in cigars):
                print >> sys.stderr, "Unexpected cigar value in ", line_number, "found", self.cigar
                print >> sys.stderr, "Found",line
                sys.exit(1)
            self.cigar_long = ""    
            for cigar in cigars:
                self.cigar_long = cigar[-1] * int(cigar[:-1])
        else:
            self.cigar_long = "*"    
            
        self.pnext_name = parts[6] 
        self.pnext_possition = int_value(line_number, line, "PNEXT", parts[7], 0, 2**31-1) 
        self.tlen = int_value(line_number, line, "TLEN", parts[8], -2**31+1, 2**31-1)
        self.seq = parts[9]
        self.qual = parts[10]

        if len(sequence_dictionary) > 0:
            sequence = sequence_dictionary[self.rname]
            if sequence:
                sequence.addAlignment(self)
            else:
                print >> sys.stderr, "Aligment line", line_number, "has an unexpected rname", self.rname
                print >> sys.stderr, "Found",line
                sys.exit(1)

default_output_path = "summary.tsv"
version = None
sorting_order = "unkown"
grouping  = "none"
               
def process_header_line(line_number, line):
    if line.startswith("@HD"):
        process_head_line(line_number, line)
    elif line.startswith("@SQ"): 
        sequence = Sequence(line_number, line)
        sequence_dictionary[sequence.name] = sequence
    elif line.startswith("@RG"): 
        print "ignoring @RG line"
    elif line.startswith("@PG"): 
        program = Program(line_number, line)
        program_dictionary[program.id] = program
    elif line.startswith("@CO"): 
        print "ignoring @PG line"
    else:
        print >> sys.stderr, "Line", line_number, "has an unexpected @ tag."
        print >> sys.stderr, "Found",line
        sys.exit(1)
        
def process_head_line(line_number, line):
    global version, sorting_order, grouping 
    if line_number == 0:
        parts = line.strip().split("\t")
        for part in parts[1:]:
            tag = part[:3]
            value = part[3:]
            if tag == "VN:":
                version = value
            elif tag == "SO:":
                if value in ["unknown","unsorted","queryname","coordinate"]:
                    sorting_order = value
                else:
                    print >> sys.stderr, "Line", line_number, "has an unexpected '@SO:' value:", value        
                    print >> sys.stderr, "Found",line
                    sys.exit(1)
            elif tag == "CO:":
                if value in ["none","query","reference"]:
                    grouping = value
                else:
                    print >> sys.stderr, "Line", line_number, "has an unexpected '@CO:' value:", value        
                    print >> sys.stderr, "Found",line
                    sys.exit(1)
            else:
                print >> sys.stderr, "Line", line_number, "has an unexpected tag",tag
                print >> sys.stderr, "Found",line
                sys.exit(1)
        if not version:        
            print >> sys.stderr, "Line", line_number, "has a no VO: tag."
            print >> sys.stderr, "Found",line
            sys.exit(1)
    else:    
        print >> sys.stderr, "Line", line_number, "has an unexpected @HD."
        print >> sys.stderr, "Found",line
        sys.exit(1)

    
def summary(input_path, 
            output_path=default_output_path):
        
    print "Summarizing ", input_path, "to", output_path
   
    in_header = True

    #read data
    with open(input_path, 'r') as input_file:
        for line_number,line in enumerate(input_file):
            if line.startswith("@"):
                if in_header:
                    process_header_line(line_number, line)
                else:
                    print >> sys.stderr, "Line", line_number, "has an unexpected @. No longer processing Header"
                    print >> sys.stderr, "Found",line
                    sys.exit(1)
            else:
                if in_header:
                    in_header == False
                Alignment(line_number, line)
                         
    with open(output_path, 'w') as output_file:
        output_file.write(strHeader())
        output_file.write("\n")
        output_file.write(str(_all_sequence))
        output_file.write("\n")
        for key in sequence_dictionary:
            sequence = sequence_dictionary[key]
            if sequence.count > 0:
                 output_file.write(str(sequence))
                 output_file.write("\n")
        
if __name__ == '__main__':
    usage = "usage: %prog [options] INPUT_FILE"
    # define options
    parser = OptionParser(usage=usage)
    parser.add_option("-o", "--output", dest="output_path", 
                      help="Path to file, where the output should be written. Defaults to counts.tsv", 
                      default="lengths.tsv")
                      
    # parse
    options, args = parser.parse_args()

    #Check exactly one input_file specified
    if len(args) < 1:
        parser.error("ERROR! No INPUT_FILE specified.")   
    elif len(args) > 1:
        parser.error("ERROR! Multiple values for INPUT_FILE found.")   
        
    summary(args[0], options.output_path)

