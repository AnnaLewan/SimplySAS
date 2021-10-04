# SimplySAS

AmpliSAS reimplementation

SimplySAS 
-----------------------------------------------------------------------------------------------------

Reimplementation of AmpliSAS. Performs genotyping of amplicon sequencing data by clustering errors and filtering artefacts following a 3-step based analysis: file preprocessing, clustering and filtering.

Files in repository:

`SimplySAS.py` - main script

`environment.yml` - environment file

`fogsaa.cpp` - alignement softwere

`primers.csv` - example file with primers data for NEW METHOD

`barcodes_primers.csv` - example file with primers and barcodes data for OLD METHOD


# To use:
1. Download git repository:

`Git clone https://github.com/AnnaLewan/Simply_SAS`

`cd SimplySAS`

2. Create environment:

`conda env create -f environment.yml`

`conda activate SimplySAS`

3. Install fogsaa

`g++ fogsaa.cpp -o fogsaa`

4. Use script as shown below:

Command:


`python SimplySAS.py -i seq_file -m new -p primer_file -o output_dir -el expected_len -se substitution_error_threshold -df min_dominant_frequency_threshold -af min_amplicon_seq_frequency`


Options:



`-i,'--INPUT',type=str,required=True, help="Input set of FASTQ or FASTA files packed into a unique .ZIP or .TAR.GZ file."`

`-m,'--METHOD',type=str,required=True, help="Method used to obtain reads (old or new). old - demultiplexing needed, new - without demultiplexing."`

`-bp, type=str, required=False, help="CSV file with primer and barcode data for old method. 5'->3' direction."`

`-p, type=str, required=False, help="CSV file with primer data for new method. 5'->3' direction."`

`-o, type=str, required=True, help="Output folder full path."`

`-el, type=int, required=False, help="Expected length of the marker sequence."`

`-se, type=float, required=False, help="Threshold for permissible substitution error (%)."`

`-df, type=float, required=False, help="Minimum frequency respect to the dominant (%)."`

`-af, type=float, required=False, help="Minimum frequency per-amplicon (%)."`

`-af, type=float, required=False, help="Minimum sequence frequency per-amplicon (%)."`

`-ch, type=int, required=False, help="Minimal length of match within sequences to consider as chimera (bp)."`

`-al, type=int, required=False, help="Maximal number of allels in one amplicon."`

`-ad, type=int, required=False, help="Minimal depth of amplicon to not be discarded."`

`-ma, type=int, required=False, help="Minimal allel depth to not be discarded in filtering."`

`-t, type=int, required=False, help="Number of threads used in analysis."`


Example:

`python SimplySAS.py -i Des_ex3.zip -m new -p 'primers.csv' -o /home/test -el 176 -se 1.56 -df 38.11 -af 43.64 -t 4`











