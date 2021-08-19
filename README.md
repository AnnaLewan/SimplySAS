# AmpliSAS_py
AmpliSAS reimplementation

PartSAS 
-----------------------------------------------------------------------------------------------------

Performs part of genotyping of amplicon sequencing data by clustering errors and filtering artefacts.

To use:
1. Download git repository:

`Git clone https://github.com/AnnaLewan/AmpliSAS_py`

`cd AmpliSAS_py`

2. Create environment:

`conda env create -f environment.yml`

`conda activate magisterka`

3. Use script as shown below: ---- NA TEN MOMENT NIE DZIAŁA (problem z implementacja klas z osobnych plików)


`python PartSAS.py -i seq_file -d csv_file -o output_dir -ml min_length -el expected_len -se substitution_error_threshold -df min_dominant_frequency_threshold -af min_amplicon_seq_frequency`

`-i,'--INPUT',type=str,required=True,help="Input set of FASTQ or FASTA files packed into a unique .ZIP or .TAR.GZ file.")`

`-d,'--DATA',type=str,required=True,help="CSV file with primer/amplicon data.")`

`-o,'--OUTPUT',type=str,required=True,help="Output folder name.")`

`-ml, '--MINLEN',type=int,required=True,help="Minimal length of sequence to consider clustering.")`

`-el, '--EXPLEN',type=int,required=True,help="Expected length of the marker sequence.")`

`-se, '--SUBERROR',type=float,required=True,help="Threshold for permissible substitution error (%).")`

`-df, '--DOMFREQ',type=float,required=True,help="Minimum frequency respect to the dominant (%).")`

`-af, '--AMFREQ',type=float,required=True,help="Minimum frequency per-amplicon (%).")`


Example:

`python PartSAS.py -i Des_ex3.zip -d startery.csv -o /test -ml 100 -el 176 -se 1.56 -df 38.11 -af 43.64`

4. Alternatywnie przez jupiter:

`jupyter notebook`

Jupyter notebook, otworzy się w przeglądarce. Nalezy wybrać plik `test.ipynb`. Parametry i ścieżki do analizowanych plików (w tym wypadku `Des_ex3.zip` i `startery.csv`) znajdują się w drugiej komórce pliku. W celu przetestowania innych danych i parametrów w drugiej komórce wpisać właściwe ścieżki i dane. W celu uruchomienia, każdą komórkę skryptu nalęży kliknąć wskazany na zdjęciu przycisk.

![PRZYCISK](/Zrzut ekranu z 2021-08-19 07-21-39.png)










