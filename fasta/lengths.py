from collections import Counter
from optparse import OptionParser, OptionGroup

class Fasta:

    def __init__(self, name):
        self.name = name
        self.sequence = []
        
    def addSequenceLine(self, line):
        self.sequence.append(line)
          
    def ___repr__(self):
        return "Fasta(" + repr(self.name) + repr(self.sequence) + ")"

    def ___str__(self):
        return ">" + self.name + "\n" + "\n".join(self.sequence) + "\n"

class FastaSummary:

    def __init__(self, fasta):
        self.name = fasta.name 
        lens = map(len,fasta.sequence)
        self.length = sum(lens)
        self.rows = len(fasta.sequence)
        self.min = min(lens)
        self.max = max(lens)
        self.counter = Counter()
        for line in fasta.sequence:
            self.counter.update(line)
 
    def ___repr__(self):
        reply = "FastaSummary(" + self.name + "," + str(self.length) + ","
        reply += str(self.rows) + "," + str(self.min) + "," + str(self.max)
        reply += repr(self.counter) + ")"
        return reply
        
    def asLine(self, keys):
        reply = str(self.length) + "\t" 
        reply += str(self.rows) + "\t" + str(self.min) + "\t" + str(self.max)
        for key in keys:            
            reply+= "\t" + str(self.counter[key]) 
        reply += "\t" + self.name 
        return reply 
        
def summary(input_path, output_path):
        
    print "Summarizing ", input_path, "to", output_path
   
    summaries = []

    #read data
    with open(input_path, 'r') as input_file:
        fasta = None
        for line in input_file:
            stripped = line.strip()
            if stripped.startswith(">"):
                if fasta:
                    summaries.append(FastaSummary(fasta))
                fasta = Fasta(stripped[1:])
            else:
                fasta.addSequenceLine(stripped)
 
    allCount = Counter()
    for summary in summaries:
        allCount.update(summary.counter)
    keys = sorted(allCount.keys())    
            
    with open(output_path, 'w') as output_file:
       line = ["Length","Rows","min Row","max_row"]
       output_file.write("\t".join(line) + "\t")
       output_file.write("\t".join(keys) + "\tName\n")
       line = [".",".",".","."]
       output_file.write("\t".join(line))
       for key in keys:            
           output_file.write("\t" + str(allCount[key])) 
       output_file.write("\tAll\n")    
       for summary in summaries:
          output_file.write(summary.asLine(keys))
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

