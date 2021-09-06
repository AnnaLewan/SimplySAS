# AmpliSAS_py
AmpliSAS reimplementation

PartSAS 
-----------------------------------------------------------------------------------------------------

Performs part of genotyping of amplicon sequencing data by clustering errors and filtering artefacts. For now it can make one cluster (from dominant sequence) and find consensus seq for this cluster.

Files in repository:

`PartSAS.py` - główny plik programu

`environment.yml` - plik umozliwiający postawienie środowiska

`fogsaa.cpp` - zewnętrzny program do alignmentu

`Salmo_MHC_I.zip` - przykładowy plik z odczytami

`test.ipynb` - plik ipythona zawierający cały skrypt i umożliwiający uruchomienie wszystkiego za pomocą jupyter notebook


# To use:
1. Download git repository:

`Git clone https://github.com/AnnaLewan/AmpliSAS_py`

`cd AmpliSAS_py`

2. Create environment:

`conda env create -f environment.yml`

`conda activate magisterka`

3. Install fogsaa

`g++ fogsaa.cpp -o fogsaa`

4. Use script as shown below:


`python PartSAS.py -i seq_file -o output_dir -ml min_length -el expected_len -se substitution_error_threshold -df min_dominant_frequency_threshold -af min_amplicon_seq_frequency`

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

`python PartSAS.py -i Des_ex3.zip -o /test -ml 100 -el 176 -se 1.56 -df 38.11 -af 43.64`

4. Alternatywnie przez jupiter:

`jupyter notebook`

Jupyter notebook, otworzy się w przeglądarce. Nalezy wybrać plik `test.ipynb`. Parametry i ścieżki do analizowanych plików znajdują się w drugiej komórce pliku. W celu przetestowania innych danych i parametrów w drugiej komórce wpisać właściwe ścieżki i dane. W celu uruchomienia, każdą komórkę skryptu nalęży kliknąć wskazany na zdjęciu przycisk.

![przycisk](/przycisk.png)










