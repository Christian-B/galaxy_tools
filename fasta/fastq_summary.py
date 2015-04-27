from collections import Counter
from optparse import OptionParser, OptionGroup

class Fastq:

    def __init__(self, name):
        self.name = name
        self.sequence = ""
        self.quality = ""
        self.plusLine = False
        
    def addLine(self, line):
        if self.plusLine:
            self.addQualityLine(line)
        else:
            self.addSequenceLine(line)
                
    def addSequenceLine(self, line):
        self.sequence+= line
          
    def addQualityLine(self, line):
        self.quality += line
        if len(self.sequence) > len(self.quality):
            print >> sys.stderr, "Unexcpected format. Quality appears longer than sequence."
            print >> sys.stderr, "Read previous fastq as:"
            print >> sys.stderr, fastq
            sys.exit(1) 
        
    def plusLineExists(self, line):
        if self.plusLine:
            return True
        self.plusLine = line
        return False   
          
    def ___repr__(self):
        reply  = "@" + self.name + "/n"
        reply += self.sequence + "/n"
        reply += "+/n"
        reply += self.quality
        return reply

_counterKeys = None

class FastqSummary:

    def __init__(self, fastq):
        if fastq:
            self.name = fastq.name 
            self.length = len(fastq.sequence)
            self.sequence_counter = Counter(fastq.sequence)
            self.minQuality = min(fastq.quality)        
            self.maxQuality = max(fastq.quality)        
            if _allSummary.name:
                _allSummary.length += self.length
                _allSummary.sequence_counter.update(self.sequence_counter)
                if _allSummary.minQuality > self.minQuality:
                    _allSummary.minQuality = self.minQuality
                if _allSummary.maxQuality < self.maxQuality:
                    _allSummary.maxQuality = self.maxQuality
                _counterKeys = None
            else:
                _allSummary.name = "All"
                _allSummary.length = self.length
                _allSummary.sequence_counter = Counter(fastq.sequence)
                _allSummary.minQuality = self.minQuality
                _allSummary.maxQuality = self.maxQuality
                _counterKeys = None
        else:
             self.name = "" 
             self.length = 0 
             self.sequence_counter = Counter()
             self.minQuality = 0
             self.maxQuality = 0
             
    def ___repr__(self):
        reply = "FastqSummary(" + self.name + "," + str(self.length) + ","
        reply += repr(self.counter) 
        reply += str(self.minQuality) + "," + str(self.maxQuality) + ")"
        return reply
        
    def __str__(self):
        reply = str(self.length) + "\t" 
        reply += str(self.minQuality) + "\t" + str(self.maxQuality)
        for key in counterKeys():            
            reply+= "\t" + str(self.sequence_counter[key]) 
        reply += "\t" + self.name 
        return reply 
      
_allSummary = FastqSummary(None)

def counterKeys():
    global _counterKeys
    if not _counterKeys:
        _counterKeys = sorted(_allSummary.sequence_counter.keys())
    return _counterKeys
            
def strHeader():
    line = ["Length","Min Q","max Q"]
    line.extend(counterKeys())
    line.append("Name")
    return "\t".join(line)
          
def summary(input_path, output_path):
        
    print "Summarizing ", input_path, "to", output_path
   
    summaries = []

    #read data
    with open(input_path, 'r') as input_file:
        fastq = None
        for (line_number, line) in enumerate(input_file):
            stripped = line.strip()
            if stripped.startswith("@"):
                if fastq:
                    if len(fastq.sequence) == len(fastq.quality):
                        #assume new Fastq
                        summaries.append(FastqSummary(fastq))
                        fastq = Fastq(stripped[1:])
                    else:
                        fastq.addQualityLine(stripped)
                else: #first fastq       
                    fastq = Fastq(stripped[1:])     
            elif stripped.startswith("+"):
                if fastq.plusLineExists(stripped):
                    fastq.addQualityLine(stripped)
            else:
                #Fastq will know if it is sequence or quality
                fastq.addLine(stripped)
        #save last Fastq       
        summaries.append(FastqSummary(fastq))
                        
    with open(output_path, 'w') as output_file:
       output_file.write(strHeader())
       output_file.write("\n")
       output_file.write(str(_allSummary))
       output_file.write("\n")
       for summary in summaries:
           output_file.write(str(summary))
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

