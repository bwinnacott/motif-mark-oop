#!/usr/bin/env python
import argparse
import cairo
import drawing_classes

GENE_POS_FACTOR = 100

def get_motifs(motif_file):
    '''
    Function to load the file containing the list of motif sequences to be searched in the 
    fasta file.
    motif_file: file containing list of motifs

    return: list of motifs
    '''
    # open connection to file
    motif = open(motif_file, 'r')
    # list for storing motif sequences
    motifs = []
    # iterate over each motif and add to list
    for line in motif:
        line = line.strip()
        motifs.append(line)
    # close the file
    motif.close()

    return motifs

def parse_fasta(fasta_file):
    '''
    Function that goes through fasta file and pulls out all relevant information for viewing 
    motifs on sequences.
    fasta_file: input fasta file containing gene sequences

    return: 
    num_seqs - number of sequences in fasta file
    longest_seq - longest sequence in fasta file
    sequences - dictionary (key: header, value: sequence)
    '''
    # initiate variables for figure dimension purposes
    num_seqs = 0
    longest_seq = 0
    # open connection to file
    with open(fasta_file,'r') as f:
        # initiate variables to store information
        header = None
        sequences = {}
        # iterate over each line
        for line in f:
            line = line.strip()
            # if header line
            if line.startswith('>'):
                num_seqs += 1
                # to account for None initially assigned to 'header'
                if header:
                    # get current longest sequence
                    longest_seq = max(longest_seq,len(sequences[header]))
                # store header
                sequences[line[1:]] = ''
                header = line[1:]
            # add the sequence line
            else:
                sequences[header] += line
        # account for the last sequence
        longest_seq = max(longest_seq,len(sequences[header]))
        
        return num_seqs,longest_seq,sequences

def get_color_palette(palette):
    '''
    Function used to get the RGB values for colors in the chosen palette.
    palette: string specified at command line for chosen color palette

    return: mapping of motif index to RGB color values
    '''
    # colorblind friendly palette
    if palette == 'colorblind_friendly':
        color_map = {0:(0.9,0.6,0),1:(0.35,0.7,0.9),2:(0,0.6,0.5),3:(0.95,0.9,0.25),4:(0,0.45,0.7)}
    # earth colors palette
    elif palette == 'earth':
        color_map = {0:(0.5,0.54,0.44),1:(0.81,0.7,0.53),2:(0.58,0.32,0.21),3:(0.16,0.27,0.25),4:(0.72,0.58,0.32)}
    # basic color palette (default)
    else:
        color_map = {0:(1,0,0),1:(0,1,0),2:(0,0,1),3:(0,1,1),4:(1,0,1)}

    return color_map
        
def output_image(input_filename,width,height,out_filetype):
    '''
    Function to write out the image to a specific file type.
    input_filename: input fasta file name
    width: width of output image
    height: height of output image
    out_filetype: file type specified at command line

    return: surface for drawing image
    '''
    # svg output
    if out_filetype == 'svg':
        return cairo.SVGSurface(str(input_filename.split('.')[0] + '.svg'),width,height)
    # png output
    elif out_filetype == 'png':
        return cairo.ImageSurface(cairo.FORMAT_ARGB32,width,height)
    # pdf output
    else:
        return cairo.PDFSurface(str(input_filename.split('.')[0] + '.pdf'),width,height)

def main():
    # set command line arguments
    parser = argparse.ArgumentParser(description='Tool used to mark motifs on gene sequences containing 1 exon and introns flanking either side. \
                                    Can handle a maximum of 5 motifs on each sequence.''')
    parser.add_argument('-f','--fa_file',action='store',required=True,type=str,help='Specifies the input FASTA file for which sequence motif marking is desired.')
    parser.add_argument('-m','--motif_file',action='store',required=True,type=str,help='Specifies the input .txt file containing the list of motifs.')
    parser.add_argument('-o','--output_type',action='store',required=True,choices=['svg','pdf','png'],type=str,help='Designates the output file type.')
    parser.add_argument('-c','--colors',action='store',choices=['earth','colorblind_friendly'],type=str,help='Specifies color palette used for differentiating motifs. \
                        If not provided at command line, the "basic" color palette is used as default.')
    # extract the argument information
    args = parser.parse_args()
    # get fasta info and motif list
    num_seq,longest_seq,sequences = parse_fasta(args.fa_file)
    motifs = get_motifs(args.motif_file)
    # get color palette and image for drawing
    color_map = get_color_palette(args.colors)
    width,height = longest_seq+10,GENE_POS_FACTOR*(num_seq+1)
    surface = output_image(args.fa_file,width,height,args.output_type)
    context = cairo.Context(surface)
    # draw the key
    key = drawing_classes.Key(motifs,color_map)
    key.draw_key(context)
    gene_ind = 1
    # draw the sequences one at a time
    for head,seq in sequences.items():
        gene = drawing_classes.Gene(seq,3)
        exon = drawing_classes.Exon(seq)
        motif_list = []
        for motif in motifs:
            motif_list.append(drawing_classes.Motif(motif,color_map,len(motifs)))
        header = drawing_classes.FastaHeader(head)
        gene_group = drawing_classes.GeneGroup(gene,exon,motif_list,header,seq,0)
        gene_group.draw(context,gene_ind)
        gene_ind += 1
    # finishing the surface object and writing out to file type specified
    if args.output_type == 'png':
        surface.write_to_png(str('Figure_1.fasta'.split('.')[0] + '.png'))
    else:
        surface.finish()

if __name__ == "__main__":
    main()