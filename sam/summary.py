from collections import Counter, OrderedDict
from optparse import OptionParser, OptionGroup
import re
import sys
       
#Based on https://samtools.github.io/hts-specs/SAMv1.pdf     
     
default_count_sequence = False
default_count_quality = False
default_count_cigar = False
default_output_path = "summary.tsv"

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

class SamHeader():

    def __init__(self,  count_sequence = default_count_sequence,      
                count_quality = default_count_quality,
                count_cigar = default_count_cigar):
        self.count_sequence = count_sequence      
        self.count_quality = count_quality
        self.count_cigar = count_cigar
        self.version = None
        self.sorting_order = "unkown"
        self.grouping  = "none"
        self.reference_dictionary = OrderedDict() 
        self.program_dictionary = OrderedDict()         
        self._all_reference = None #To be added later
        self._needs_updating= True
     
    def all_reference(self):
        if not self._all_reference:
            self._all_reference = Reference(None, "All", self)
        if not "*" in self.reference_dictionary:                   
            self.reference_dictionary["*"] = Reference(None, "Unkown", self)        
        return self._all_reference
                    
    def update(self):
        all_reference = self.all_reference()
        if self.count_sequence:    
            self._sequence_keys = sorted(all_reference.sequence_counter.keys())
        if self.count_quality:    
            self._quality_keys = sorted(all_reference.quality_counter.keys())
        if self.count_cigar:
            self._cigar_keys = sorted(all_reference.cigar_counter.keys())
        self._needs_updating = False
        
    def sequenceKeys(self):
        if self._needs_updating:
            self.update()
        return self._sequence_keys
            
    def qualityKeys(self):
        if self._needs_updating:
            self.update()
        return self._quality_keys
            
    def cigarKeys(self):
        if self._needs_updating:
            self.update()
        return self._cigar_keys
        
    def header_line(self):
        line = ["Count","Min Q","max Q"]
        if self.count_sequence:    
            line.extend(self.sequenceKeys())
            line.append("TSeq")
        if self.count_quality:
            line.extend(self.qualityKeys())
            line.append("TQual")
        if self.count_cigar:
            line.extend(self.cigarKeys())    
            line.append("Tcigar")
        line.append("Name")
        return "\t".join(line)  
     
    def addAlignment(self, alignment):
        allreference = self.all_reference()
        if len(self.reference_dictionary) > 0:
            reference = self.reference_dictionary[alignment.rname]
            if reference:
                reference.addAlignment(alignment)
            else:
                print >> sys.stderr, "Aligment line", line_number, "has an unexpected rname", alignment.rname
                print >> sys.stderr, "Found",line
                sys.exit(1)
        allreference.addAlignment(alignment)
        self._needs_updating= True
     
    def process_header_line(self, line_number, line):
        if line.startswith("@HD"):
            self.process_head_line(line_number, line)
        elif line.startswith("@SQ"): 
            reference = Reference(line_number, line, self)
            self.reference_dictionary[reference.name] = reference
        elif line.startswith("@RG"): 
            print "ignoring Read Group(@RG) line"
        elif line.startswith("@PG"): 
            program = Program(line_number, line)
            self.program_dictionary[program.id] = program
        elif line.startswith("@CO"): 
            print "ignoring Comment(@CO) line"
        else:
            print >> sys.stderr, "Line", line_number, "has an unexpected @ tag."
            print >> sys.stderr, "Found",line
            sys.exit(1)
        
    def process_head_line(self, line_number, line):
        if line_number == 0:
            parts = line.strip().split("\t")
            for part in parts[1:]:
                tag = part[:3]
                value = part[3:]
                if tag == "VN:":
                    self.version = value
                elif tag == "SO:":
                    if value in ["unknown","unsorted","queryname","coordinate"]:
                        self.sorting_order = value
                    else:
                        print >> sys.stderr, "Line", line_number, "has an unexpected '@SO:' value:", value        
                        print >> sys.stderr, "Found",line
                        sys.exit(1)
                elif tag == "CO:":
                    if value in ["none","query","reference"]:
                        self.grouping = value
                    else:
                        print >> sys.stderr, "Line", line_number, "has an unexpected '@CO:' value:", value        
                        print >> sys.stderr, "Found",line
                        sys.exit(1)
                else:
                    print >> sys.stderr, "Line", line_number, "has an unexpected tag",tag
                    print >> sys.stderr, "Found",line
                    sys.exit(1)
            if not self.version:        
                print >> sys.stderr, "Line", line_number, "has a no VO: tag."
                print >> sys.stderr, "Found",line
                sys.exit(1)
        else:    
            print >> sys.stderr, "Line", line_number, "has an unexpected @HD."
            print >> sys.stderr, "Found",line
            sys.exit(1)
            
    def __str__(self):
        reply = "count_sequence: " + str(self.count_sequence)       
        reply+= "count_quality: " + str(self.count_quality)
        reply+= "count_cigar: " + str(self.count_cigar)
        return reply 
                         
class Reference:
    def __init__(self, line_number, line, header=SamHeader()):
        self.header = header
        if self.header.count_sequence:
            self.sequence_counter = Counter()
        if self.header.count_quality:    
            self.quality_counter = Counter()
        if self.header.count_cigar:   
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
        if self.header.count_sequence:
            self.sequence_counter.update(alignment.seq)
        if self.header.count_quality:  
            self.quality_counter.update(alignment.qual)
        if self.header.count_cigar:   
            self.cigar_counter.update(alignment.cigar_long)
        if self.max_mapq:
            if self.max_mapq < alignment.mapq:
                self.max_mapq = alignment.mapq       
            if self.min_mapq > alignment.mapq:
                self.min_mapq = alignment.mapq       
        else:
            self.max_mapq = alignment.mapq        
            self.min_mapq = alignment.mapq       
 
    def __str__(self):
        reply = str(self.count) 
        reply+=  "\t" 
        reply+= str(self.min_mapq) 
        reply+= "\t" 
        reply+= str(self.max_mapq)
        if self.header.count_sequence:
            for key in self.header.sequenceKeys():            
                reply+= "\t"
                reply+= str(self.sequence_counter[key]) 
            reply+= "\t" 
            reply+= str(sum(self.sequence_counter.values()) )  
        if self.header.count_quality:  
            for key in self.header.qualityKeys():            
                reply+= "\t"
                reply+= str(self.quality_counter[key]) 
            reply+= "\t" 
            reply+= str(sum(self.quality_counter.values()))   
        if self.header.count_cigar:  
            for key in self.header.cigarKeys():            
                reply+= "\t"
                reply+= str(self.cigar_counter[key]) 
            reply+= "\t" 
            reply+= str(sum(self.cigar_counter.values()))   
        reply+= "\t" 
        reply+= self.name 
        return reply      

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
    def __init__(self, line_number, line, header = SamHeader()):
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
        
        self.optional_fields = {}
        for optional in parts[11:]:
            tag = optional[:2]
            the_type = optional[3]
            value = optional[5:]
            if tag in self.optional_fields:
                print >> sys.stderr, "Found duplicate optional tag", tag, " in line",line_number
                print >> sys.stderr, "Found",line
                sys.exit(1)
            if the_type == "A": #Printable character
                if len(value_string) != 1:
                    print >> sys.stderr, "Found duplicate optional tag", tag, " in line",line_number
                    print >> sys.stderr, "Found",line
                    sys.exit(1)
                else:
                    self.optional_fields[tag] = value
            elif the_type == "i": #Signed integer
                try:
                    self.optional_fields[tag] = int(value)
                except:
                    print >> sys.stderr, "Line", line_number, "has an none ingeter", tag, "value"
                    print >> sys.stderr, "Found",line
                    sys.exit(1)
            elif the_type == "f": #Single-precision Floating number
                try:
                    self.optional_fields[tag] = float(value)
                except:
                    print >> sys.stderr, "Line", line_number, "has an none float", tag, "value"
                    print >> sys.stderr, "Found",line
                    sys.exit(1)
            elif the_type == "Z": #Printable string, including space
                self.optional_fields[tag] = value
            elif the_type == "H": #Byte array in the Hex format
                try:
                    self.optional_fields[tag] = bytearray.frmohex(value)
                except:
                    print >> sys.stderr, "Line", line_number, "has an none hex/bytearray", tag, "value"
                    print >> sys.stderr, "Found",line
                    sys.exit(1)
            elif the_type == "B": #Integer or numeric array
                #todo parse string
                self.optional_fields[tag] = value    
            else:    
                print >> sys.stderr, "Line", line_number, "has an unexpected type", the_type
                print >> sys.stderr, "Found",line
                sys.exit(1)

        #Now that we have the key information add the alignment to the summary
        header.addAlignment(self)
    
def summary(input_path, 
            output_path=default_output_path,
            count_sequence=default_count_sequence,
            count_quality=default_count_quality,
            count_cigar=default_count_cigar):

    print "Summarizing ", input_path, "to", output_path

    in_header = True
    header = SamHeader(count_sequence, count_quality, count_cigar)
    
    print header
    
    #read data
    with open(input_path, 'r') as input_file:
        for line_number,line in enumerate(input_file):
            if line.startswith("@"):
                if in_header:
                    header.process_header_line(line_number, line)
                else:
                    print >> sys.stderr, "Line", line_number, "has an unexpected @. No longer processing Header"
                    print >> sys.stderr, "Found",line
                    sys.exit(1)
            else:
                if in_header:
                    in_header == False
                Alignment(line_number, line, header)
                         
    with open(output_path, 'w') as output_file:
        output_file.write(header.header_line())
        output_file.write("\n")
        output_file.write(str(header.all_reference()))
        output_file.write("\n")
        for key in header.reference_dictionary:
            reference = header.reference_dictionary[key]
            if reference.count > 0:
                 output_file.write(str(reference))
                 output_file.write("\n")
        
if __name__ == '__main__':
    usage = "usage: %prog [options] INPUT_FILE"
    # define options
    parser = OptionParser(usage=usage)
    parser.add_option("-o", "--output", dest="output_path", 
                      help="Path to file, where the output should be written. Defaults to counts.tsv", 
                      default="lengths.tsv")
    count_group = OptionGroup(parser, "Count individual Values Options",
                    "Specifes if the number of times each values is found in these fields. "
                    "For example it will count how many of each of C, G, T, and A are found in the sequence. "
                    "Does not effect validation of lengths.")
    count_group.add_option("-s", "--sequence", 
            dest="count_sequence",
            action="store_true",
            help="If used this flag will count the values in every aligmnent's sequence",
            default = default_count_sequence)
    count_group.add_option("-q", "--quality", 
            dest="count_quality",
            action="store_true",
            help="If used this flag will count the values in every aligmnent's quality",
            default = default_count_quality)
    count_group.add_option("-c", "--cigar", 
            dest="count_cigar",
            action="store_true",
            help="If used this flag will count the values in every aligmnent's sequence",
            default = default_count_cigar)
    parser.add_option_group(count_group)
                      
    # parse
    options, args = parser.parse_args()

    #Check exactly one input_file specified
    if len(args) < 1:
        parser.error("ERROR! No INPUT_FILE specified.")   
    elif len(args) > 1:
        parser.error("ERROR! Multiple values for INPUT_FILE found.")   
        
    summary(args[0], options.output_path, options.count_sequence,
            options.count_quality, options.count_cigar)

