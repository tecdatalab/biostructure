import csv
import os


def write_in_file(file_path, headers, data):
    real_file_path = os.path.abspath(file_path)
    if os.path.isfile(real_file_path):
        with open(real_file_path, 'a', encoding="utf-8") as outfile:
            writer = csv.writer(outfile)
            for datum in data:
                writer.writerow(datum)
    else:
        with open(file_path, 'w', encoding="utf-8") as outfile:
            writer = csv.writer(outfile)
            writer.writerow(headers)
            for datum in data:
                writer.writerow(datum)
