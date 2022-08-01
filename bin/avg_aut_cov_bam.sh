../bin/samtools-1.11/samtools idxstats ${1} | sed s/"chr"/""/g | awk '$1<=38 {print $0}' | grep -v "*" | awk '$2>0{print $3/$2}' | awk 'BEGIN{s=0;}{s=s+$1;}END{print s/NR;}'
