"""Compress and move every file in folder `_generated_data` to the versioned folders `annotations` and `genome_info`"""

import os
import zipfile

# Function to read file contents
def read_file(file_path):
    with open(file_path, 'rb') as f:
        return f.read()

# Function to read a file inside a zip archive
def read_file_from_zip(zip_path, file_name_in_zip):
    with zipfile.ZipFile(zip_path, 'r') as zipf:
        with zipf.open(file_name_in_zip) as f:
            return f.read()

for dir_path, _, filenames in os.walk('_generated_data'):

    if len(filenames) == 0:
        continue

    destination_dir = dir_path.replace('_generated_data', '.', 1)

    os.makedirs(destination_dir, exist_ok=True)

    for filename in filenames:
        source_file = dir_path  + '/' + filename
        dest_file = destination_dir  + '/' + filename + '.zip'

        if os.path.exists(dest_file):
            source_file_content = read_file(source_file)
            dest_file_content = read_file_from_zip(dest_file, filename)
            if source_file_content == dest_file_content:
                continue

        print(f'creating {dest_file} from {source_file}')
        with zipfile.ZipFile(dest_file, 'w', zipfile.ZIP_DEFLATED) as zipf:
            zipf.write(source_file, arcname=filename)




