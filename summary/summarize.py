from optparse import OptionParser, OptionGroup
import sys

default_output_path = "summary.tsv"
default_seperator = "\t"
default_header = "-1" #No header line
default_data_start = 0 #Top row (header line always ignored)

def summerize(input_path, 
              output_path=default_output_path, 
              seperator=default_seperator, 
              header_line=default_header, 
              data_start=default_data_start):
           
    def processDataLine(parts):    
        for (col, st) in enumerate(parts): 
            if st:
                counts[col]+= 1
                if not types[col]:
                    types[col] = "int"
                    st_maxs[col] = st
                    st_mins[col] = st
                    first = True
                else:    
                    first = False
                    if (st > st_maxs[col]):
                        st_maxs[col] = st
                    if (st < st_mins[col]):
                        st_mins[col] = st
                if types[col] == "int":
                    try:
                        value = int(st)
                    except ValueError:
                        types[col] = "float"
                if types[col] == "float":
                    try:
                        value= float(st)
                    except ValueError:
                        types[col] = "string"
                if types[col] != "string":
                    sums[col]+= value
                    if first or (value > maxs[col]):
                         maxs[col] = value
                    if first or (value < mins[col]):
                         mins[col] = value
                     
    print "Summarizing ", input_path, "to", output_path

    if header_line >= 0:
        first_line = min(header_line, data_start)
    else:
        first_line = data_start    
        
    with open(input_path, 'r') as input_file:
        for line_number,line in enumerate(input_file):
            if line_number >= data_start or line_number == header_line:
                parts = line.strip().split(seperator) 
                if line_number == first_line:    
                    number_of_columns = len(parts)
                    counts = [0 for i in range(number_of_columns)]
                    sums = [0 for i in range(number_of_columns)]
                    st_maxs = [None for i in range(number_of_columns)]
                    st_mins = [None for i in range(number_of_columns)]
                    maxs = [None for i in range(number_of_columns)]
                    mins = [None for i in range(number_of_columns)]
                    types = [None for i in range(number_of_columns)]
                    headers = []
                elif len(parts) != number_of_columns: 
                    if line_number >= header_line and line_number >= data_start:
                        print >> sys.stderr, "Line", line_number, "does not have the expected number of columns."
                        print >> sys.stderr, "Found",parts,"but expected",number_of_columns,"column(s)."
                        sys.exit(1)
                    #else save to ignore as not a header or data row
                if line_number == header_line:
                    headers.extend(parts)   
                else:
                    processDataLine(parts)

    with open(output_path, 'w') as output_file:
        line = ["Column","Type","count","Min","Max","Sum","Average"]
        output_file.write("\t".join(line) + "\n")
        for col in range(number_of_columns):    
            line = []
            if headers:
                line.append(headers[col])
            else:
                line.append("column " + str(col) + ":")
            line.append(types[col])
            line.append(str(counts[col]))
            if types[col] == "string":
                line.append(st_mins[col])
                line.append(st_maxs[col])
                line.append("n/a")
                line.append("n/a")
                line.append("")
            else:
                line.append(str(mins[col]))
                line.append(str(maxs[col]))
                line.append(str(sums[col]))
                line.append(str(sums[col]/float(counts[col])))
            output_file.write("\t".join(line) + "\n")
        
if __name__ == '__main__':
    usage = "usage: %prog [options] INPUT_FILE"
    # define options
    parser = OptionParser(usage=usage)
    parser.add_option("-o", "--output", dest="output_path", 
                      help="Path to file, where the output should be written. Defaults to " + default_output_path, 
                      default=default_output_path)
    line_group = OptionGroup(parser, "Line Options",
                    "Specifes the Lines numbers of Header and Data. "
                    "All line numbers are 0 based so the top line is 0. "
                    "If specified a header line can be below the first data line.")
    line_group.add_option("--header", dest="header_line",
                      help="Line number of Header Line. Default is (-1) no header.",
                      default = default_header)
    line_group.add_option("-d", "--data_start", dest="data_start",
                      help="Line number of First line of Data. Default is " + str(default_data_start), 
                      default = default_data_start)
    parser.add_option_group(line_group)
    seperator_group = OptionGroup(parser, "Seperator Options",
                    "Specifes the Seperator to used to split columns. "
                    "Defaults to tab. "
                    "Only one Seperator option can be provided.")
    seperator_group.add_option("-a", "--ascii_seperator", dest="ascii_seperator", help="ASCii code of seperator")
    seperator_group.add_option("-s", "--seperator", dest="seperator", help="Column seperator as character.")
    parser.add_option_group(seperator_group)

    # parse
    options, args = parser.parse_args()

    #Check exactly one input_file specified
    if len(args) < 1:
        parser.error("ERROR! No INPUT_FILE specified.")   
    elif len(args) > 1:
        parser.error("ERROR! Multiple values for INPUT_FILE found.")   
 
    #Get and check the column seperator
    if options.ascii_seperator:
        if options.seperator:
            parser.error("options ascii_seperator and seperator are mutually exclusive")
        else:
            seperator = chr(int(options.ascii_seperator))  
    else:
        if options.seperator:
            seperator = options.seperator
        else:
            seperator = default_seperator
            
    #Check line options
    try:
        header_line = int(options.header_line)
    except ValueError:
        parser.error("ERROR! None Integer header: " + options.header_line)   
    try:
        data_start = int(options.data_start)
        if data_start < 0:
            parser.error("ERROR! data_start must be 0 or positive.")   
    except ValueError:
        parser.error("ERROR! None Integer data_start: " + options.data_start)   
        
    
    summerize(args[0], options.output_path, seperator, header_line, data_start)

