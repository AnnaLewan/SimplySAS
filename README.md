# SimplySAS

AmpliSAS reimplementation

SimplySAS 
-----------------------------------------------------------------------------------------------------

Reimplementation of AmpliSAS. Performs genotyping of amplicon sequencing data by clustering errors and filtering artefacts following a 3-step based analysis: de-multiplexing, clustering and filtering.

Files in repository:

`SimplySAS.py` - main script

`environment.yml` - environment file

`fogsaa.cpp` - alignement softwere

`Salmo_MHC_I.zip` - example multifile for NEW METHOD

`primers.csv` - example file with primers data for NEW METHOD

`barcodes_primers.csv` - example file with primers and barcodes data for OLD METHOD


# To use:
1. Download git repository:

`Git clone https://github.com/AnnaLewan/AmpliSAS_py`

`cd AmpliSAS_py`

2. Create environment:

`conda env create -f environment.yml`

`conda activate SimplySAS`

3. Install fogsaa

`g++ fogsaa.cpp -o fogsaa`

4. Use script as shown below:


`python PartSAS.py -i seq_file -o output_dir -el expected_len -se substitution_error_threshold -df min_dominant_frequency_threshold -af min_amplicon_seq_frequency`

`-i,'--INPUT',type=str,required=True,help="Input set of FASTQ or FASTA files packed into a unique .ZIP or .TAR.GZ file."`

`-o,'--OUTPUT',type=str,required=True,help="Output folder name."`

`-el, '--EXPLEN',type=int,required=False,help="Expected length of the marker sequence."`

`-se, '--SUBERROR',type=float,required=False,help="Threshold for permissible substitution error (%)."`

`-df, '--DOMFREQ',type=float,required=False,help="Minimum frequency respect to the dominant (%)."`

`-af, '--AMFREQ',type=float,required=False,help="Minimum frequency per-amplicon (%)."`

`-af", '--AMFREQ',type=float,required=False,help="Minimum sequence frequency per-amplicon (%)."`

`-ch", '--CHIMERA',type=int,required=False,help="Minimal length of match within sequences to consider as chimera (bp)."`

`-al", '--MAXALL',type=int,required=False,help="Maximal number of allels in one amplicon."`

`-ad", '--MINDEPTH',type=int,required=False,help="Minimal depth of amplicon to not be discarded."`

`-mf", '--MINFREQ',type=float,required=False,help="Minimal frequency of consensus sequence (allel) to not be discarded in filtering."`

`-nc", '--NONCOD',type=int,required=False,help="Discard noncoding sequences (1 -on; default off)."`


Example:

`python PartSAS.py -i Des_ex3.zip -o /test -el 176 -se 1.56 -df 38.11 -af 43.64`

4. Alternatywnie przez jupiter:

`jupyter notebook`

Jupyter notebook, otworzy się w przeglądarce. Nalezy wybrać plik `test.ipynb`. Parametry i ścieżki do analizowanych plików znajdują się w drugiej komórce pliku. W celu przetestowania innych danych i parametrów w drugiej komórce wpisać właściwe ścieżki i dane. W celu uruchomienia, każdą komórkę skryptu nalęży kliknąć wskazany na zdjęciu przycisk.

![przycisk](/przycisk.png)










