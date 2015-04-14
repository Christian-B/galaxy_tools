from optparse import OptionParser, OptionGroup
import sys

def summerize(input_path, output_path, seperator, has_header):
        
    def processLine(line):
        for (col, st) in enumerate(line.strip().split(seperator)): 
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
                        print value
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

    with open(input_path, 'r') as input_file:
        header = input_file.readline().strip()
        parts = header.split(seperator) 
        number_of_columns = len(parts)
        counts = [0 for i in range(number_of_columns)]
        sums = [0 for i in range(number_of_columns)]
        st_maxs = [None for i in range(number_of_columns)]
        st_mins = [None for i in range(number_of_columns)]
        maxs = [None for i in range(number_of_columns)]
        mins = [None for i in range(number_of_columns)]
        types = [None for i in range(number_of_columns)]
        if has_header:
            headers = parts
        else:
            processLine(header)
        for line in input_file:
            processLine(line)

    with open(output_path, 'w') as output_file:
        line = ["Column","Type","count","Min","Max","Sum","Average"]
        output_file.write("\t".join(line) + "\n")
        for col in range(number_of_columns):    
            line = []
            if has_header:
                line.append(headers[col])
            else:
                line.append("column " + str(col) + ":")
            line.append(types[col])
            line.append(str(counts[col]))
            if types[col] == "string":
                line.append(st_mins[col])
                line.append(st_maxs[col])
                line.append("not numerical")
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
                      help="Path to file, where the output should be written. Defaults to summary.tsv", 
                      default="summary.tsv")
    parser.add_option("--header", action="store_true", dest="has_header",
                      help="Considers the first line a header line. Default is no header.")
    group = OptionGroup(parser, "Seperator Options",
                    "Specifes the Seperator to used to split columns. "
                    "Defaults to tab. "
                    "Only one Seperator option can be provided.")
    group.add_option("-a", "--ascii_seperator", dest="ascii_seperator", help="ASCii code of seperator")
    group.add_option("-s", "--seperator", dest="seperator", help="Column seperator as character.")
    parser.add_option_group(group)

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
            parser.error("options -a and -b are mutually exclusive")
        else:
            seperator = chr(int(options.ascii_seperator))  
    else:
        if options.seperator:
            seperator = options.seperator
        else:
            seperator = "\t"
    
    summerize(args[0], options.output_path, seperator, options.has_header)

