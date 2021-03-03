#!/usr/bin/env python
import cairo

LEFT_MARGIN = 5
GENE_POS_FACTOR = 100
GENE_POS_OFFSET = 50
EXON_POS_OFFSET = 25

def check_degen_symb(symbol):
    '''
    Function that takes in an IUPAC base symbol and returns the base(s) which it represents.
    symbol: IUPAC base symbol

    return: list of matching nts for input symbol
    '''
    # set the symbol mapping
    bases = {
        'A':['A','a'],'a':['A','a'],
        'C':['C','c'],'c':['C','c'],
        'G':['G','g'],'g':['G','g'],
        'T':['T','U','t','u'],'t':['T','U','t','u'],
        'U':['U','T','u','t'],'u':['U','T','u','t'],
        'W':['A','T','a','t'],'w':['A','T','a','t'],
        'S':['C','G','c','g'],'s':['C','G','c','g'],
        'M':['A','C','a','c'],'m':['A','C','a','c'],
        'K':['G','T','g','t'],'k':['G','T','g','t'],
        'R':['A','G','a','g'],'r':['A','G','a','g'],
        'Y':['C','T','c','t'],'y':['C','T','c','t'],
        'B':['C','G','T','c','g','t'],'b':['C','G','T','c','g','t'],
        'D':['A','G','T','a','g','t'],'d':['A','G','T','a','g','t'],
        'H':['A','C','T','a','c','t'],'h':['A','C','T','a','c','t'],
        'V':['A','C','G','a','c','g'],'v':['A','C','G','a','c','g'],
        'N':['A','C','G','T','a','c','g','t'],'n':['A','C','G','T','a','c','g','t'],
        'Z':[''],'z':['']
        }

    return bases[symbol]

def chars_to_skip(motif):
    '''
    This function takes in a motif pattern and generates an array containing the number of 
    characters to skip for each index position in the pattern, to optimize matching. The array 
    is returned. This is the preprocessing step for the KMP pattern matching algorithm.
    motif: sequence motif

    return: list of characters to skip for each letter in motif (when matching new window)
    '''
    # get length of motif
    l = len(motif)
    # initiate list of 0's of size l
    skip_chars = [0] * l
    # set counter to track index of current character in 'motif'; set another to track length 
    # of previous longest prefix that is also a suffix in 'motif' (i.e., matching sub-patterns)
    # Note: we start the first counter at 1, due to there not being a proper prefix for a single character
    i = 1
    prev_len = 0
    # loop until all characters are accounted for in 'motif'
    while i < l:
        # if characters match, add 1 to length counter and assign the current index in 'skip_chars' to 
        # the new length; increment the current index
        if motif[i] == motif[prev_len]:
            prev_len += 1
            skip_chars[i] = prev_len
            i += 1
        else:
            # if the longest prefix length is currently 0, assign current index in 'skip_chars' to 0, 
            # increment index counter
            if prev_len == 0:
                skip_chars[i] = 0
                i += 1
            # if previous longest length of prefix is not 0, assign new length (i.e., new pattern index for 
            # which to start comparing new window)
            else:
                prev_len = skip_chars[prev_len-1]

    return skip_chars

def search_sequence(motif,sequence):
    '''
    Function using KMP pattern matching algorithm to search a sequence for all instances of 
    current motif. All positions for which the motif is found are returned.
    motif: sequence motif
    sequence: sequence in fasta file

    return: list of indexes where a given motif matches the in the gene sequence
    '''
    # initiate list to store motif mapping positions
    motif_pos = []
    # get lengths of motif and sequence
    len_pat = len(motif)
    len_seq = len(sequence)
    # call 'chars_to_skip' function to get count of characters to be skipped
    skip_chars = chars_to_skip(motif)
    # initiate counters for sequence and motif index
    i = 0
    j = 0
    # loop until all characters are evaluated in sequence
    while i < len_seq:
        # if characters match, increment counters
        if sequence[i] in check_degen_symb(motif[j]):
            i += 1
            j += 1
        # if motif matches current window in sequence, get matching starting position; get 
        # index of next character in motif to be matched
        if j == len_pat:
            motif_pos.append(i-j)
            j = skip_chars[j-1]
        # if current characters don't match
        elif i < len_seq and sequence[i] not in check_degen_symb(motif[j]):
            # if first character of motif doesn't match window, shift window by 1
            if j == 0:
                i += 1
            # if current character (other than first) doesn't match, get index for next 
            # character in motif to be matched
            else:
                j = skip_chars[j-1]

    return motif_pos

class Key:

    def __init__(self,motifs,color_map):
        self.motifs = motifs
        self.color_map = color_map
        self.curr_key_len = 200
        self.col_ind = 0

    def draw_motif(self,context,motif):
        context.rectangle(self.curr_key_len,25,15,30)
        color = self.color_map[self.col_ind]
        (r,g,b) = color
        context.set_source_rgb(r,g,b)
        context.fill()
        # adjust position of text
        self.curr_key_len += 25
        context.move_to(self.curr_key_len,45)
        context.set_source_rgb(0,0,0)
        # add motif sequence
        context.show_text(motif)
        # depending on length of motif, adjust position of next motif label
        if len(motif) <= 2:
            self.curr_key_len += 25
        elif 2<len(motif)<=4:
            self.curr_key_len += 40
        elif 4<len(motif)<=6:
            self.curr_key_len += 55
        elif 6<len(motif)<=8:
            self.curr_key_len += 70
        else:
            self.curr_key_len += 85

    def draw_key(self,context):
        context.rectangle(15,25,30,30)
        context.fill()
        context.move_to(65,45)
        context.show_text('Exon')
        context.set_line_width(3)
        context.move_to(105,42)
        context.line_to(135,42)
        context.stroke()
        context.move_to(155,45)
        context.show_text('Intron')
        for m in self.motifs:
            self.draw_motif(context,m)
            self.col_ind += 1
        context.rectangle(LEFT_MARGIN,10,self.curr_key_len-5,60)
        context.stroke()

class Gene:

    def __init__(self,sequence,line_width):
        self.length = len(sequence)
        self.line_width = line_width

    def draw(self,context,gene_ind):
        # draw the intron on the surface
        context.set_line_width(self.line_width)
        context.move_to(LEFT_MARGIN,gene_ind*GENE_POS_FACTOR+GENE_POS_OFFSET)
        context.line_to(self.length+LEFT_MARGIN,gene_ind*GENE_POS_FACTOR+GENE_POS_OFFSET)
        context.set_source_rgb(0,0,0)
        context.stroke()

class Exon:

    def __init__(self,sequence):
        self.sequence = [i for i,char in enumerate(sequence) if char.isupper()]
        self.start_pos = min(self.sequence)
        self.end_pos = max(self.sequence)

    def draw(self,context,gene_ind):
        # draw the exon on the surface
        context.rectangle(self.start_pos+LEFT_MARGIN,EXON_POS_OFFSET+GENE_POS_FACTOR*gene_ind,self.end_pos-self.start_pos+1,GENE_POS_OFFSET)
        context.set_source_rgb(0,0,0)
        context.fill()

class Motif:

    def __init__(self,motif,color_map,num_motifs):
        self.length = len(motif)
        self.color_map = color_map
        self.num_motifs = num_motifs
        self.motif = motif

    def draw(self,context,gene_ind,motif_ind,motif_map):
        color = self.color_map[motif_ind]
        (r,g,b) = color
        # get factors for adjusting the position of the motif on the image
        pos_factor = 0 if self.num_motifs == 1 else 25/(self.num_motifs-1)
        len_factor = 50 if self.num_motifs == 1 else 25
        # draw the motif for each location it matched to the gene sequence
        for pos in motif_map[self.motif]:
            context.rectangle(pos+LEFT_MARGIN,EXON_POS_OFFSET+(pos_factor*motif_ind)+GENE_POS_FACTOR*gene_ind,self.length,len_factor)
            context.set_source_rgba(r,g,b)
            context.fill()

class FastaHeader:

    def __init__(self,text):
        self.text = text

    def draw(self,context,gene_ind):
        # set position on surface for header
        context.set_source_rgb(0,0,0)
        context.move_to(LEFT_MARGIN,gene_ind*GENE_POS_FACTOR+10)
        # modify the header to make it more readable
        header = self.text.split(' ')
        context.show_text('Gene: ' + str(header[0]) + 2*' ' + '-->' + 2*' ' + 
                        'Chromosome: ' + header[1].split(':')[0] + 2*' ' + 
                        '-->' + 2*' ' + 'Coordinates: ' + str(header[1].split(':')[1]))

class GeneGroup:

    def __init__(self,Gene,Exon,Motif_list,FastaHeader,seq,motif_ind):
        self.Gene = Gene
        self.Exon = Exon
        self.Motif_list = Motif_list
        self.FastaHeader = FastaHeader
        self.motif_ind = motif_ind
        # initiate dictionary
        self.motif_mapping = {}
        # loop over each motif and add motif as key and list of indexes in sequence where 
        # it matches as value
        for motif in self.Motif_list:
            self.motif_mapping[motif.motif] = search_sequence(motif.motif,seq)

    def draw(self,context,gene_ind):
        self.Gene.draw(context,gene_ind)
        self.Exon.draw(context,gene_ind)
        for motif in self.Motif_list:
            motif.draw(context,gene_ind,self.motif_ind,self.motif_mapping)
            self.motif_ind += 1
        self.FastaHeader.draw(context,gene_ind)