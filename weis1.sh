#!/bin/bash

i=1 n=0
while IFS=, read line; do
    col1=$(echo ${line} | cut -d , -f 1 | tr -d "\r\n" | sed -e 's/^"//' -e 's/"$//')
    ((n >= i)) && \
    mkdir /tmp/fastq_files/arr_seq/${col1} && \
    mkdir /tmp/fastq_files/arr_seq/${col1}/${col1}_1 && \
    mkdir /tmp/fastq_files/arr_seq/${col1}/${col1}_2
    ((n++))
done < weis_metadata_1210.csv && \
while IFS=, read line; do
    col1=$(echo ${line} | cut -d , -f 1 | tr -d "\r\n" | sed -e 's/^"//' -e 's/"$//')
    col2=$(echo ${line} | cut -d , -f 2 | tr -d "\r\n" | sed -e 's/^"//' -e 's/"$//')
    gzip /tmp/fastq_files/ena_files/${col1}_1.fastq.gz --decompress --stdout > /tmp/fastq_files/arr_seq/${col2}/${col2}_1/${col1}_1.fastq && \
    gzip /tmp/fastq_files/ena_files/${col1}_2.fastq.gz --decompress --stdout > /tmp/fastq_files/arr_seq/${col2}/${col2}_2/${col1}_2.fastq
done < id_list.csv && \
while IFS=, read line; do
    col2=$(echo ${line} | cut -d , -f 2 | tr -d "\r\n" | sed -e 's/^"//' -e 's/"$//')
    list1=$(eval "ls /tmp/fastq_files/arr_seq/${col2}/${col2}_1/*.fastq")
    list2=$(eval "ls /tmp/fastq_files/arr_seq/${col2}/${col2}_2/*.fastq")
    for i in ${list1}; do cat $i >> /tmp/fastq_files/arr_seq/${col2}/${col2}_1.fastq ; done
    for i in ${list2}; do cat $i >> /tmp/fastq_files/arr_seq/${col2}/${col2}_2.fastq ; done
done < id_list.csv && \
i=1 n=0
while IFS=, read line; do
    col1=$(echo ${line} | cut -d , -f 1 | tr -d "\r\n" | sed -e 's/^"//' -e 's/"$//')
    ((n >= i)) && \
    gzip /tmp/fastq_files/arr_seq/${col1}/${col1}_1.fastq && \
    gzip /tmp/fastq_files/arr_seq/${col1}/${col1}_2.fastq
    ((n++))
done < weis_metadata_1210.csv && \
i=1 n=0
while IFS=, read line; do
    col1=$(echo ${line} | cut -d , -f 1 | tr -d "\r\n" | sed -e 's/^"//' -e 's/"$//')
    ((n >= i)) && \
    rm -r /tmp/fastq_files/arr_seq/${col1}/${col1}_1 && \
    rm -r /tmp/fastq_files/arr_seq/${col1}/${col1}_2
    ((n++))
done < weis_metadata_1210.csv && \
echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" > manifest.tsv && \
i=1 n=0
while IFS=, read line; do
    col1=$(echo ${line} | cut -d , -f 1 | tr -d "\r\n" | sed -e 's/^"//' -e 's/"$//')
    ((n >= i)) && \
    echo -e "${col1}\t"'$PWD'"/fastq_files/arr_seq/${col1}/${col1}_1.fastq.gz\t"'$PWD'"/fastq_files/arr_seq/${col1}/${col1}_2.fastq.gz" >> manifest.tsv
    ((n++))
done < weis_metadata_1210.csv
