from optparse import OptionParser, OptionGroup
import sys

def countFrequency(input_path, output_path, seperator, max_options, has_header):
        
    def processLine(line):
        for (col, st) in enumerate(line.strip().split(seperator)): 
            if len(counts[col]) <= max_options:
                  counts[col][st] = counts[col].setdefault(st, 0) + 1                      

    print "Counting ", input_path, "to", output_path

    with open(input_path, 'r') as input_file:
        header = input_file.readline().strip()
        parts = header.split(seperator) 
        number_of_columns = len(parts)
        counts = [{} for i in range(number_of_columns)]
        if has_header:
            headers = parts
        else:
            processLine(header)
        for line in input_file:
            processLine(line)

    with open(output_path, 'w') as output_file:
        line = ["Column","Value","Count"]
        output_file.write("\t".join(line) + "\n")
        for col in range(number_of_columns):
            if len(counts[col]) <= max_options:   
                for value in sorted(counts[col]):
                    line = []
                    if has_header:
                        line.append(headers[col])
                    else:
                        line.append("column " + str(col) + ":")
                    line.append(value)
                    line.append(str(counts[col][value]))
                    output_file.write("\t".join(line) + "\n")
        
if __name__ == '__main__':
    usage = "usage: %prog [options] INPUT_FILE"
    # define options
    parser = OptionParser(usage=usage)
    parser.add_option("-o", "--output", dest="output_path", 
                      help="Path to file, where the output should be written. Defaults to counts.tsv", 
                      default="counts.tsv")
    parser.add_option("--header", action="store_true", dest="has_header",
                      help="Considers the first line a header line. Default is no header.")
    parser.add_option("-m", "--maxoptions", dest="max_options", default="10",
                      help="Maximum number of different values to count per column. Must be a positive integer. Default is 10.")
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
    
    #Get and check max Options
    try:
        max_options = int(options.max_options)
        if max_options <= 0:
            parser.error("ERROR! maxoptions must be greater than 0.")   
    except ValueError:
        parser.error("ERROR! None Integer maxoptions: " + options.max_options)   

    
    
    countFrequency(args[0], options.output_path, seperator, max_options, options.has_header)

