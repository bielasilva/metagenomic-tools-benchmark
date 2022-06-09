#!/bin/bash

while getopts 'f:d:ph' flag; do
    case "${flag}" in
        f) FILES="${OPTARG}" ;;
        d) SAVE_DIR="${OPTARG}" ;;
        p) PROTEIN='true';;
        * | h) printf 'USAGE: genome_downloader.sh -f ASSEMBLY_SUMMARY -d SAVE_DIR (-p)\n' &&  exit 1;;
    esac
done

SCRIPT_DIR=$(dirname "$0")
MODE="$MODE"

if [[ $PROTEIN ]]; then
    MODE="*gbff.gz"
fi

[ ! -d $SAVE_DIR ] && mkdir $SAVE_DIR 

echo "Iniciando análise"
cat $FILES > assembly_tot.txt

i=1

TOTAL=$(awk -F "\t" '$12=="Complete Genome" && $11=="latest"{print $20}' assembly_tot.txt | wc -l)
echo "Total para download: $TOTAL "

until [ `ls $SAVE_DIR/$MODE | wc -l` == $TOTAL ]; do
    echo "Tentativa $i"
    echo "Checando já baixados"
    ls $SAVE_DIR/$MODE > downed.txt
    awk -F "\t" '$12=="Complete Genome" && $11=="latest"' assembly_tot.txt > filtered.txt
    $SCRIPT_DIR/check_downloaded.py filtered.txt downed.txt > checked.tmp
    
    if [[ $PROTEIN ]]; then
        awk 'BEGIN{FS="\t";OFS="/"}$12=="Complete Genome" && $11=="latest"{l=split($20,a,"/");print $20,a[l]"_genomic.gbff.gz"}' checked.tmp > list_down.txt
    else
        awk 'BEGIN{FS="\t";OFS="/"}$12=="Complete Genome" && $11=="latest"{l=split($20,a,"/");print $20,a[l]"_genomic.fna.gz"}' checked.tmp > list_down.txt
    fi
    
    FALTA=$(cat list_down.txt | wc -l)
    echo "Faltam $FALTA"
    
    echo "Tentando normal"
    cat list_down.txt | xargs -P 3 -n 1 wget -c -P $SAVE_DIR -nv 2> err.tmp
    cat list_down.txt | xargs -P 3 -n 1 wget -c -P $SAVE_DIR -nv 2> err.tmp
    
    if [[ `ls $SAVE_DIR/$MODE | wc -l` == $TOTAL ]]; then
        break
    fi
    
    echo "Tentando inverso"
    tac list_down.txt | xargs -P 3 -n 1 wget -c -P $SAVE_DIR -nv 2> err.tmp
    tac list_down.txt | xargs -P 3 -n 1 wget -c -P $SAVE_DIR -nv 2> err.tmp
    echo ""

    ((i=i+1))
    if [ $i == "150" ]; then
        echo "Download não terminado"
        break
    fi
done

if [ `ls $SAVE_DIR/$MODE | wc -l` == $TOTAL ]; then
    echo "Download terminado com sucesso"
    read -p "Deseja remover os arquivos intermediários? y/n  " choice
    if [ $choice == "y" ]; then
        rm downed.txt checked.tmp err.tmp assembly_tot.txt list_down.txt filtered.txt
    elif [ $choice == "n" ]; then
        echo "OK, aquivos mantidos"
    fi
fi
