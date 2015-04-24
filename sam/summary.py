from collections import Counter, OrderedDict
from optparse import OptionParser, OptionGroup
import re
import sys
       
#Based on https://samtools.github.io/hts-specs/SAMv1.pdf     
     
default_level = "ignore"
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

    def __init__(self):
        self.version = None
        self.sorting_order = "unkown"
        self.grouping  = "none"
        self.reference_dictionary = OrderedDict() 
        self.program_dictionary = OrderedDict()         
        self._all_reference = None #To be added later
        self._needs_updating= True
     
    def all_reference(self):
        if not self._all_reference:
            self._all_reference = Reference(None, "All")
        if not "*" in self.reference_dictionary:                   
            self.reference_dictionary["*"] = Reference(None, "Unkown")        
        return self._all_reference
                    
    def update(self):
        all_reference = self.all_reference()
        self._flag_keys = sorted(all_reference.flag_counter.keys())
        self._cigar_keys = sorted(all_reference.cigar_counter.keys())
        self._mapq_keys = sorted(all_reference.mapq_counter.keys())
        self._sequence_keys = sorted(all_reference.sequence_counter.keys())
        self._quality_keys = sorted(all_reference.quality_counter.keys())
        self._needs_updating = False
        
    def flag_keys(self):
        if self._needs_updating:
            self.update()
        return self._flag_keys

    def cigar_keys(self):
        if self._needs_updating:
            self.update()
        return self._cigar_keys
        
    def mapq_keys(self):
        if self._needs_updating:
            self.update()
        return self._mapq_keys
        
    def sequence_keys(self):
        if self._needs_updating:
            self.update()
        return self._sequence_keys
            
    def quality_keys(self):
        if self._needs_updating:
            self.update()
        return self._quality_keys
            
    def addAlignment(self, alignment, line_number, line):
        allreference = self.all_reference()
        if len(self.reference_dictionary) > 0:
            if alignment.rname in self.reference_dictionary:
                reference = self.reference_dictionary[alignment.rname]
                reference.addAlignment(alignment)
            else:
                print self.reference_dictionary.keys()
                print >> sys.stderr, "Aligment line", line_number, "has an unexpected rname", alignment.rname
                print >> sys.stderr, "Found",line
                sys.exit(1)
        allreference.addAlignment(alignment)
        self._needs_updating= True
     
    def process_header_line(self, line_number, line):
        if line.startswith("@HD"):
            self.process_head_line(line_number, line)
        elif line.startswith("@SQ"): 
            reference = Reference(line_number, line)
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
                                     
class Reference:
    def __init__(self, line_number, line):
        self.flag_counter = Counter()
        self.count_position = 0       
        self.sum_position = 0                    
        self.max_position = 0        
        self.min_position = 0  
        self.mapq_counter = Counter()
        self.cigar_count = 0                    
        self.cigar_counter = Counter()
        self.sequence_counter = Counter()
        self.quality_counter = Counter()
        self.length = None
        self.count = 0
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
        self.flag_counter[alignment.flag]+=1
        if alignment.position > 0:
            if self.max_position:
                if self.max_position < alignment.position:
                    self.max_position = alignment.position       
                if self.min_position > alignment.position:
                    self.min_positions = alignment.position       
                self.count_position+= 1        
                self.sum_position+= alignment.position                    
            else:
                self.count_position = 1        
                self.sum_position = alignment.position                   
                self.max_position = alignment.position        
                self.min_position = alignment.position                     
        self.mapq_counter[alignment.mapq]+=1
        if alignment.cigar_long != "*":
            self.cigar_count+= 1
            self.cigar_counter.update(alignment.cigar_long)
        self.sequence_counter.update(alignment.sequence)
        self.quality_counter.update(alignment.quality)
 
    @staticmethod
    def header_line(flag=False, position=False, mapq=False, cigar=False, 
                    sequence=False, quality=False):
        line = ["Count"]
        if flag:
            if flag == True:
                line.extend(["FlagMin","FlagMax"])
            else:
                for key in flag:
                    line.append(str(key))
        if position:
            line.extend(["#Pos","PosMin","PosMax","posAvg"])
        if mapq:
            line.append("#mapQ")
            if mapq == True:
                line.extend(["mapQMin","mapQMax","mapQAvg"])
            else:
                for key in mapq:
                    line.append(str(key))
        if cigar:
            line.append("#cigar")
            if cigar != True:
                line.extend(cigar)    
            line.append("LenCigar")
            line.append("AvgLen")
        if sequence:
            if sequence != True:
                line.extend(sequence)    
            line.append("LenSeq")
        if quality:
            if quality == True:
                line.extend(["QualMin","QualMax"])
            else:
                line.extend(quality)    
            line.append("LenQual")            
        line.append("Name")
        return "\t".join(line)  
     
    def summary_line(self, flag=False, position=False, mapq=False, cigar=False, 
                     sequence=False, quality=False):
        line = [str(self.count)]
        if flag:
            if flag == True:
                line.append(str(min(self.flag_counter)))
                line.append(str(max(self.flag_counter)))
            else:
                for key in flag:
                    line.append(str(self.flag_counter[key]))
        if position:
            line.append(str(self.count_position))   
            line.append(str(self.max_position))
            line.append(str(self.min_position))                      
            if self.count_position == 0:
                line.append("0")
            else:        
                line.append(str(self.sum_position/float(self.count_position)))
        if mapq:
            count = sum(self.mapq_counter.values())
            line.append(str(count))
            if mapq == True:
                line.append(str(min(self.mapq_counter)))
                line.append(str(max(self.mapq_counter)))
                total = sum(map(lambda x: x[0] * x[1], self.mapq_counter.items()))
                line.append(str(total/float(count)))
            else:
                for key in mapq:
                    line.append(str(self.mapq_counter[key]))
        if cigar:
            line.append(str(self.cigar_count))
            if cigar != True:
                for key in cigar:
                    line.append(str(self.cigar_counter[key]))
            length = sum(self.cigar_counter.values())       
            line.append(str(length))
            if self.cigar_count > 0:
                line.append(str(length/float(self.cigar_count)))
            else:    
                line.append("0")
        if sequence:
            if sequence != True:
                for key in sequence:
                    line.append(str(self.sequence_counter[key]))
            line.append(str(sum(self.sequence_counter.values())))
        if quality:
            if quality == True:
                line.append(str(min(self.quality_counter)))
                line.append(str(max(self.quality_counter)))
            else:    
                for key in quality:
                    line.append(str(self.quality_counter[key]))
            line.append(str(sum(self.quality_counter.values())))
        line.append(self.name)
        return "\t".join(line)  

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
      
        self.flag = int_value(line_number, line, "FLAG", parts[1], 0, 2**16-1)
        self.multiple_segments = self.flag & 1 != 0
        self.properly_aligned = self.flag & 2 != 0
        self.unmapped = self.flag & 4 != 0
        self.next_unmapped = self.flag & 4 != 0
        self.reverse_complemented = self.flag & 16 != 0  #0x10
        self.next_reverse_complemented = self.flag & 32 != 0 #OX20
        self.first_segment = self.flag & 64 != 0 #Ox40
        self.last_segment = self.flag & 128 != 0 #Ox80
        self.secondary = self.flag & 256 != 0 #Ox100
        self.failed_quality_control = self.flag & 512 != 0 #Ox200
        self.duplicate = self.flag & 1024 != 0 #Ox400
        self.supplement = self.flag & 2048 != 0 #Ox800
      
        self.rname = parts[2]              
        self.position = int_value(line_number, line, "POS", parts[3], 0, 2**31-1)
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
        self.sequence = parts[9]
        self.quality = parts[10]
        
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
    
def processOption(parameter, name, keys):
    if parameter.lower() == "ignore":
        return None
    elif parameter.lower() == "summary":
        return True   
    elif parameter.lower() == "count":
        if keys:
            return keys
        else:    
            print >> sys.stderr, "Count value not supperted for",name,"parameter"
    else:
        print >> sys.stderr, "Unexpected value for",name,"parameter. Found", parameter
        sys.exit(1)

def summary(input_path, 
            output_path=default_output_path,
            flag=default_level,
            position=default_level,
            mapq=default_level,
            cigar=default_level,
            sequence=default_level,
            quality=default_level,
            flag_path=None,
            mapq_path=None,
            cigar_path=None,
            sequence_path=None,
            quality_path=None):

    print "Summarizing ", input_path, "to", output_path

    in_header = True
    header = SamHeader()
        
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
                alignment = Alignment(line_number, line)
                header.addAlignment(alignment, line_number, line)
    
    flag_show = processOption(flag,"flag", header.flag_keys())
    position_show = processOption(position,"position", None)
    mapq_show = processOption(mapq,"mapq", header.mapq_keys())
    cigar_show = processOption(cigar,"mapq", header.cigar_keys())
    sequence_show = processOption(sequence,"sequence", header.sequence_keys())
    quality_show = processOption(quality,"quality", header.quality_keys())
        
    with open(output_path, 'w') as output_file:
        output_file.write(Reference.header_line(flag=flag_show,
                                                position=position_show,
                                                mapq=mapq_show,
                                                cigar=cigar_show,
                                                sequence=sequence_show,
                                                quality=quality_show))
        output_file.write("\n")
        output_file.write(str(header.all_reference().summary_line(flag=flag_show,
                                                position=position_show,
                                                mapq=mapq_show,
                                                cigar=cigar_show,
                                                sequence=sequence_show,
                                                quality=quality_show)))
        output_file.write("\n")
        for key in header.reference_dictionary:
            reference = header.reference_dictionary[key]
            if reference.count > 0:
                 output_file.write(reference.summary_line(flag=flag_show,
                                                position=position_show,
                                                mapq=mapq_show,
                                                cigar=cigar_show,
                                                sequence=sequence_show,
                                                quality=quality_show))
                 output_file.write("\n")
                 
    if flag_path:               
        with open(flag_path, 'w') as output_file:
            output_file.write(Reference.header_line(flag=header.flag_keys()))
            output_file.write("\n")
            output_file.write(str(header.all_reference().summary_line(flag=header.flag_keys())))
            output_file.write("\n")
            for key in header.reference_dictionary:
                reference = header.reference_dictionary[key]
                if reference.count > 0:
                     output_file.write(reference.summary_line(flag=header.flag_keys()))
                     output_file.write("\n")
    if mapq_path:               
        with open(mapq_path, 'w') as output_file:
            output_file.write(Reference.header_line(mapq=header.mapq_keys()))
            output_file.write("\n")
            output_file.write(str(header.all_reference().summary_line(mapq=header.mapq_keys())))
            output_file.write("\n")
            for key in header.reference_dictionary:
                reference = header.reference_dictionary[key]
                if reference.count > 0:
                     output_file.write(reference.summary_line(mapq=header.mapq_keys()))
                     output_file.write("\n")
    if cigar_path:               
        with open(cigar_path, 'w') as output_file:
            output_file.write(Reference.header_line(cigar=header.cigar_keys()))
            output_file.write("\n")
            output_file.write(str(header.all_reference().summary_line(cigar=header.cigar_keys())))
            output_file.write("\n")
            for key in header.reference_dictionary:
                reference = header.reference_dictionary[key]
                if reference.count > 0:
                     output_file.write(reference.summary_line(cigar=header.cigar_keys()))
                     output_file.write("\n")

    if sequence_path:               
        with open(sequence_path, 'w') as output_file:
            output_file.write(Reference.header_line(sequence=header.sequence_keys()))
            output_file.write("\n")
            output_file.write(str(header.all_reference().summary_line(sequence=header.sequence_keys())))
            output_file.write("\n")
            for key in header.reference_dictionary:
                reference = header.reference_dictionary[key]
                if reference.count > 0:
                     output_file.write(reference.summary_line(sequence=header.sequence_keys()))
                     output_file.write("\n")
                     
    if quality_path:               
        with open(quality_path, 'w') as output_file:
            output_file.write(Reference.header_line(quality=header.quality_keys()))
            output_file.write("\n")
            output_file.write(str(header.all_reference().summary_line(quality=header.quality_keys())))
            output_file.write("\n")
            for key in header.reference_dictionary:
                reference = header.reference_dictionary[key]
                if reference.count > 0:
                     output_file.write(reference.summary_line(quality=header.quality_keys()))
                     output_file.write("\n")
        
if __name__ == '__main__':
    usage = "usage: %prog [options] INPUT_FILE"
    # define options
    parser = OptionParser(usage=usage)
    parser.add_option("-o", "--output", dest="output_path", 
                      help="Path to general file, where any output not redirected below should be written. Defaults to summary.tsv", 
                      default="summary.tsv")
    data_group = OptionGroup(parser, "Summary/Count individual Values Options",
                    "Specifes what action should be taken with the data for this column. "
                    "Legal values are 'ignore'(default), 'summary', 'count'."
                    "'ignore' or not provided will cause the column in this data not to be output."
                    "'summary' will cause summary columns to be added to the output file."
                    "'count' will cause the counts of each value found to be added to the output file."
                    "For example it will count how many of each of C, G, T, and A are found in the sequence. "
                    "Does not effect validation of lengths.")
    data_group.add_option("-f", "--flag", 
            dest="flag",
            help="Causes the flag column to be output as requested",
            default = default_level)
    data_group.add_option("-p", "--position", 
            dest="position",
            help="Causes the position column to be output as requested. No Count position available",
            default = default_level)
    data_group.add_option("-m", "--mapq", 
            dest="mapq",
            help="Causes the mapQ column to be output as requested",
            default = default_level)
    data_group.add_option("-c", "--cigar", 
            dest="cigar",
            help="Causes the cigar column to be output as requested",
            default = default_level)
    data_group.add_option("-s", "--sequence", 
            dest="sequence",
            help="Causes the sequence column to be output as requested",
            default = default_level)
    data_group.add_option("-q", "--quality", 
            dest="quality",
            help="Causes the quality column to be output as requested",
            default = default_level)
    parser.add_option_group(data_group)
    
    path_group = OptionGroup(parser, "Path to files to write Count details to"
           "If provided will create a second output file with the counts of the values found this column"
           "Not effected by the Summary/Count flags."
           "Must be a file path and not a directory path")       
    path_group.add_option("--flag_path", 
            dest="flag_path",
            help="Causes the counts of flag column to be output in this path",
            default = None)
    path_group.add_option("--mapq_path", 
            dest="mapq_path",
            help="Causes the counts of mapQ column to be output in this path",
            default = None)
    path_group.add_option("--cigar_path", 
            dest="cigar_path",
            help="Causes the counts of cigar column to be output in this path",
            default = None)
    path_group.add_option("--sequence_path", 
            dest="sequence_path",
            help="Causes the counts of sequence column to be output in this path",
            default = None)
    path_group.add_option("--quality_path", 
            dest="quality_path",
            help="Causes the counts of quality column to be output in this path",
            default = None)
    parser.add_option_group(path_group)
                      
    # parse
    options, args = parser.parse_args()

    #Check exactly one input_file specified
    if len(args) < 1:
        parser.error("ERROR! No INPUT_FILE specified.")   
    elif len(args) > 1:
        parser.error("ERROR! Multiple values for INPUT_FILE found.")   
        
    summary(args[0], 
           options.output_path, 
           flag=options.flag, 
           position=options.position,
           mapq=options.mapq,
           cigar=options.cigar, 
           sequence=options.sequence,
           quality=options.quality,
           flag_path=options.flag_path,
           mapq_path=options.mapq_path,
           cigar_path=options.cigar_path,
           sequence_path=options.sequence_path,
           quality_path=options.quality_path)

