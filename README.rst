Welcome to illuminator-data
===========================

You'll find here the scripts used to generate the data required by illuminator, and the data itself. The data is generated from R package SeSAMe.
The manifest is directly created from data downloaded `here <https://zwdzwd.github.io/InfiniumAnnotation>`_, and the genome information is extracted from the R objects.

The scripts ``create_manifest.py`` and ``create_genome_infos.r`` will generate the data as .csv files in the folder ``generated_data``.
The versioned and compressed data used by illuminator is found in folders ``annotations``, ``arrays``, and ``genome_info``.
To move the generated data to the versioned folders and compress it, run the ``update_data.py`` script.

